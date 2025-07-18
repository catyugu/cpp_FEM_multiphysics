// src/core/coupling/ElectroThermalCoupling.cpp
#include "core/coupling/ElectroThermalCoupling.hpp"
#include "core/mesh/Element.hpp"
#include "core/mesh/TriElement.hpp"
#include "core/mesh/TetElement.hpp"
#include "core/sources/VolumetricSource.hpp"
#include "physics/PhysicsField.hpp"
#include "physics/Current1D.hpp"
#include "physics/Current2D.hpp"
#include "physics/Current3D.hpp"
#include "utils/SimpleLogger.hpp"
#include <string>
#include <core/mesh/LineElement.hpp>
#include <core/FEValues.hpp> // Include FEValues

namespace Core {

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

        const std::string joule_heat_tag = "joule_heating_source";
        heat_field_->removeSourcesByTag(joule_heat_tag);

        auto& logger = Utils::Logger::instance();
        logger.info("    Applying Joule Heat as a tagged volumetric source...");

        const auto* mesh = emag_field_->getMesh();
        const auto* dof_manager = emag_field_->getDofManager();
        const auto& material = emag_field_->getMaterial();
        const auto& emag_solution = emag_field_->getSolution();
        const auto& heat_solution = heat_field_->getSolution();

        // Determine the dimensionality of the problem from the type of the physics field.
        // This is a robust way to ensure we only process the correct element types.
        bool is_3d = dynamic_cast<Physics::Current3D*>(emag_field_) != nullptr;
        bool is_2d = dynamic_cast<Physics::Current2D*>(emag_field_) != nullptr;
        bool is_1d = dynamic_cast<Physics::Current1D*>(emag_field_) != nullptr;

        for (const auto& elem_ptr : mesh->getElements()) {
            double T_avg = 293.15; // Default to room temp
            if (heat_solution.size() > 0) {
                // FIX: Get all DOFs for the heat field on this element, respecting its order
                // Ensure element's mathematical order is correctly applied before getting DOFs
                elem_ptr->setOrder(heat_field_->getElementOrder());
                const auto heat_element_dofs = heat_field_->get_element_dofs(elem_ptr);

                T_avg = 0.0;
                size_t valid_dofs_count = 0;
                for (int dof_idx : heat_element_dofs) {
                    if (dof_idx != -1) { // Check if DOF exists
                        T_avg += heat_solution(dof_idx);
                        valid_dofs_count++;
                    }
                }
                if (valid_dofs_count > 0) {
                    T_avg /= valid_dofs_count;
                } else {
                    T_avg = 293.15; // Fallback if no valid DOFs found (shouldn't happen with correct setup)
                }
            }

            const double local_sigma = material.getProperty("electrical_conductivity", T_avg);
            double total_element_power = 0.0;

            // --- Dimension-Specific Calculations ---
            if (is_3d) {
                if (auto* tet_elem = dynamic_cast<TetElement*>(elem_ptr)) {
                    tet_elem->setOrder(emag_field_->getElementOrder()); // Ensure element has correct order
                    auto fe_values = tet_elem->create_fe_values(emag_field_->getElementOrder()); // Use emag_field's order for quadrature as well for consistency

                    double integrated_joule_heat_density = 0.0;
                    for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                        fe_values->reinit(q_p);
                        const auto& B_matrix_from_fe = fe_values->get_shape_gradients();
                        const double detJ_x_w = fe_values->get_detJ_times_weight();

                        // Get element DOFs for the EMag field, respecting its order
                        const auto element_dofs = emag_field_->get_element_dofs(tet_elem);
                        Eigen::VectorXd nodal_voltages(tet_elem->getNumNodes());
                        for (size_t k = 0; k < tet_elem->getNumNodes(); ++k) {
                            if (element_dofs[k] != -1) {
                                nodal_voltages(k) = emag_solution(element_dofs[k]);
                            } else {
                                nodal_voltages(k) = 0.0;
                            }
                        }

                        Eigen::Vector3d grad_V = B_matrix_from_fe * nodal_voltages;
                        integrated_joule_heat_density += local_sigma * grad_V.squaredNorm() * detJ_x_w;
                    }
                    total_element_power = integrated_joule_heat_density;
                }
            } else if (is_2d) {
                if (auto* tri_elem = dynamic_cast<TriElement*>(elem_ptr)) {
                    tri_elem->setOrder(emag_field_->getElementOrder());
                    auto fe_values = tri_elem->create_fe_values(emag_field_->getElementOrder());

                    double integrated_joule_heat_density = 0.0;
                    for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                        fe_values->reinit(q_p);
                        const auto& B_matrix_from_fe = fe_values->get_shape_gradients();
                        const double detJ_x_w = fe_values->get_detJ_times_weight();

                        const auto element_dofs = emag_field_->get_element_dofs(tri_elem);
                        Eigen::VectorXd nodal_voltages(tri_elem->getNumNodes());
                        for (size_t k = 0; k < tri_elem->getNumNodes(); ++k) {
                            if (element_dofs[k] != -1) {
                                nodal_voltages(k) = emag_solution(element_dofs[k]);
                            } else {
                                nodal_voltages(k) = 0.0;
                            }
                        }

                        Eigen::Vector2d grad_V = B_matrix_from_fe * nodal_voltages;
                        integrated_joule_heat_density += local_sigma * grad_V.squaredNorm() * detJ_x_w;
                    }
                    total_element_power = integrated_joule_heat_density;
                }
            } else if (is_1d) {
                if (auto* line_elem = dynamic_cast<LineElement*>(elem_ptr)) {
                    line_elem->setOrder(emag_field_->getElementOrder());
                    auto fe_values = line_elem->create_fe_values(emag_field_->getElementOrder());

                    double integrated_joule_heat_density = 0.0;
                    for (size_t q_p = 0; q_p < fe_values->num_quadrature_points(); ++q_p) {
                        fe_values->reinit(q_p);
                        const auto& B_matrix_from_fe = fe_values->get_shape_gradients();
                        const double detJ_x_w = fe_values->get_detJ_times_weight();

                        const auto element_dofs = emag_field_->get_element_dofs(line_elem);
                        Eigen::VectorXd nodal_voltages(line_elem->getNumNodes());
                        for (size_t k = 0; k < line_elem->getNumNodes(); ++k) {
                            if (element_dofs[k] != -1) {
                                nodal_voltages(k) = emag_solution(element_dofs[k]);
                            } else {
                                nodal_voltages(k) = 0.0;
                            }
                        }
                        Eigen::Vector<double, 1> grad_V = B_matrix_from_fe * nodal_voltages;
                        integrated_joule_heat_density += local_sigma * grad_V.squaredNorm() * detJ_x_w;
                    }
                    total_element_power = integrated_joule_heat_density;
                }
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