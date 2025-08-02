#include "core/coupling/ElectroThermalCoupling.hpp"
#include "core/mesh/Element.hpp"
#include "physics/PhysicsField.hpp"
#include "utils/SimpleLogger.hpp"
#include <string>
#include <core/FEValues.hpp>
#include "core/ReferenceElement.hpp"

namespace Core {
    // ... (setup function remains the same) ...
    void ElectroThermalCoupling::setup(std::vector<Physics::PhysicsField*>& fields) {
        auto& logger = Utils::Logger::instance();
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

        const auto* mesh = emag_field_->getMesh();
        const auto& material = emag_field_->getMaterial();
        const auto& emag_solution = emag_field_->getSolution();
        const auto& heat_solution = heat_field_->getSolution();

        Eigen::VectorXd& F_coupling = heat_field_->getCouplingRHS();
        if (F_coupling.size() != heat_field_->getRHS().size()) {
            F_coupling.resize(heat_field_->getRHS().size());
        }
        F_coupling.setZero();

        for (const auto& elem_ptr : mesh->getElements()) {
            double total_element_power = 0.0;
            int quad_order = 2;

            const auto& emag_ref_data = Core::ReferenceElementCache::get(elem_ptr->getTypeName(), elem_ptr->getNodes().size(), emag_field_->getElementOrder(), quad_order);
            Core::FEValues emag_fe_values(elem_ptr->getGeometry(), emag_field_->getElementOrder(), emag_ref_data);

            const auto& heat_ref_data = Core::ReferenceElementCache::get(elem_ptr->getTypeName(), elem_ptr->getNodes().size(), heat_field_->getElementOrder(), quad_order);
            Core::FEValues heat_fe_values(elem_ptr->getGeometry(), heat_field_->getElementOrder(), heat_ref_data);

            const auto emag_dofs = emag_field_->getElementDofs(elem_ptr);
            const auto heat_dofs = heat_field_->getElementDofs(elem_ptr);

            Eigen::VectorXd nodal_voltages(emag_dofs.size());
            for (size_t k = 0; k < emag_dofs.size(); ++k) nodal_voltages(k) = (emag_dofs[k] != -1) ? emag_solution(emag_dofs[k]) : 0.0;

            Eigen::VectorXd nodal_temperatures(heat_dofs.size());
            for (size_t k = 0; k < heat_dofs.size(); ++k) nodal_temperatures(k) = (heat_dofs[k] != -1) ? heat_solution(heat_dofs[k]) : 0.0;

            double integrated_joule_heat = 0.0;
            for (size_t q_p = 0; q_p < emag_fe_values.num_quadrature_points(); ++q_p) {
                emag_fe_values.reinit(q_p);
                heat_fe_values.reinit(q_p);

                const auto& emag_B = emag_fe_values.get_shape_gradients();
                const auto& heat_N = heat_fe_values.get_shape_values();
                const double detJ_x_w = emag_fe_values.get_detJ_times_weight();

                double temp_at_qp = heat_N.transpose() * nodal_temperatures;
                const double sigma_at_qp = material.getProperty("electrical_conductivity", {{"Temperature", temp_at_qp}});

                Eigen::VectorXd grad_V = emag_B * nodal_voltages;
                integrated_joule_heat += sigma_at_qp * grad_V.squaredNorm() * detJ_x_w;
            }
            total_element_power = integrated_joule_heat;

            if (total_element_power > 0) {
                std::vector<double> integral_Ni_dV(heat_dofs.size(), 0.0);
                double total_element_volume = 0.0;
                for (size_t q_p = 0; q_p < heat_fe_values.num_quadrature_points(); ++q_p) {
                    heat_fe_values.reinit(q_p);
                    const auto& N = heat_fe_values.get_shape_values();
                    const double detJ_x_w = heat_fe_values.get_detJ_times_weight();
                    for (size_t i = 0; i < N.size(); ++i) integral_Ni_dV[i] += N(i) * detJ_x_w;
                    total_element_volume += detJ_x_w;
                }

                if (std::abs(total_element_volume) > 1e-12) {
                    for (size_t i = 0; i < heat_dofs.size(); ++i) {
                        if (heat_dofs[i] != -1) {
                            F_coupling(heat_dofs[i]) += total_element_power * (integral_Ni_dV[i] / total_element_volume);
                        }
                    }
                }
            }
        }
    }

} // namespace Core