#include "core/coupling/ElectroThermalCoupling.hpp"
#include "core/mesh/Element.hpp"
#include "core/mesh/TriElement.hpp"
#include "core/sources/VolumetricSource.hpp"
#include "physics/PhysicsField.hpp"
#include "utils/SimpleLogger.hpp"
#include <string>

namespace Core {

    void ElectroThermalCoupling::setup(std::vector<Physics::PhysicsField*>& fields) {
        auto& logger = SimpleLogger::Logger::instance();
        logger.info("Setting up Electro-Thermal coupling...");

        for (auto* field : fields) {
            if (std::string(field->getVariableName()) == "Voltage") {
                emag_field_ = field;
            } else if (std::string(field->getVariableName()) == "Temperature") {
                heat_field_ = field;
            }
        }

        if (!emag_field_ || !heat_field_) {
            logger.error("Could not find both EMag and Heat fields for coupling.");
        } else {
            logger.info("Successfully coupled EMag and Heat fields.");
        }
    }

    void ElectroThermalCoupling::execute() {
        if (!emag_field_ || !heat_field_) return;

        const std::string joule_heat_tag = "joule_heating_source";

        // Step 1: Safely remove all sources from the previous iteration using the tag.
        heat_field_->removeSourcesByTag(joule_heat_tag);

        auto& logger = SimpleLogger::Logger::instance();
        logger.info("    Applying Joule Heat as a tagged volumetric source...");

        const auto& mesh = emag_field_->getMesh();
        const auto& dof_manager = emag_field_->getDofManager();
        const auto& material = emag_field_->getMaterial();
        const auto& emag_solution = emag_field_->getSolution();
        const auto& heat_solution = heat_field_->getSolution();

        // Step 2: Calculate and add new sources with the specific tag.
        for (const auto& elem_ptr : mesh->getElements()) {
            double T_avg = 300.0;
            if (heat_solution.size() > 0) {
                T_avg = 0.0;
                for (const auto& node : elem_ptr->getNodes()) {
                    int dof_idx = dof_manager->getEquationIndex(node->getId(), "Temperature");
                    if (dof_idx != -1) {
                        T_avg += heat_solution(dof_idx);
                    }
                }
                T_avg /= elem_ptr->getNumNodes();
            }

            const double local_sigma = material.getProperty("electrical_conductivity", T_avg);
            double total_element_power = 0.0;

            // Calculate the total power (Watts) generated in the element
            if (auto* tri_elem = dynamic_cast<TriElement*>(elem_ptr)) {
                auto B = tri_elem->getBMatrix();
                auto nodes = tri_elem->getNodes();
                Eigen::Vector3d nodal_voltages;
                for (int j = 0; j < 3; ++j) {
                    nodal_voltages(j) = emag_solution(dof_manager->getEquationIndex(nodes[j]->getId(), "Voltage"));
                }
                Eigen::Vector2d grad_V = B * nodal_voltages;
                double joule_heat_density = local_sigma * grad_V.squaredNorm(); // Power density (W/m^3)
                double element_volume = tri_elem->getArea() * 1.0; // Assume 1.0m thickness for 2D
                total_element_power = joule_heat_density * element_volume;
            } else if (auto* line_elem = dynamic_cast<LineElement*>(elem_ptr)) {
                auto nodes = line_elem->getNodes();
                double h = line_elem->getLength();
                double V_i = emag_solution(dof_manager->getEquationIndex(nodes[0]->getId(), "Voltage"));
                double V_j = emag_solution(dof_manager->getEquationIndex(nodes[1]->getId(), "Voltage"));
                double E = std::abs(V_i - V_j) / h;
                double joule_heat_density = local_sigma * E * E;
                double element_volume = h * 1.0;
                total_element_power = joule_heat_density * element_volume;
            }

            // Add the source term to the heat field if power is generated
            if (total_element_power > 0) {
                auto source = std::make_unique<VolumetricSource>(
                    elem_ptr->getId(), total_element_power, joule_heat_tag
                );
                heat_field_->addSource(std::move(source));
            }
        }
    }

} // namespace Core