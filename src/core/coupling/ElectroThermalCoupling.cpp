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
#include <core/FEValues.hpp> // 包含 FEValues

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
        const auto& material = emag_field_->getMaterial();
        const auto& emag_solution = emag_field_->getSolution();
        const auto& heat_solution = heat_field_->getSolution();

        bool is_3d = dynamic_cast<Physics::Current3D*>(emag_field_) != nullptr;
        bool is_2d = dynamic_cast<Physics::Current2D*>(emag_field_) != nullptr;
        bool is_1d = dynamic_cast<Physics::Current1D*>(emag_field_) != nullptr;

        for (const auto& elem_ptr : mesh->getElements()) {
            double total_element_power = 0.0;

            // --- 分维度进行精确计算 ---
            if (is_3d) {
                if (auto* tet_elem = dynamic_cast<TetElement*>(elem_ptr)) {
                    // 为电场和热场分别创建FEValues对象，以处理它们可能不同的单元阶次
                    tet_elem->setOrder(emag_field_->getElementOrder());
                    auto emag_fe_values = tet_elem->create_fe_values(emag_field_->getElementOrder());

                    tet_elem->setOrder(heat_field_->getElementOrder());
                    auto heat_fe_values = tet_elem->create_fe_values(heat_field_->getElementOrder());

                    const auto emag_dofs = emag_field_->getElementDofs(tet_elem);
                    const auto heat_dofs = heat_field_->getElementDofs(tet_elem);

                    Eigen::VectorXd nodal_voltages(emag_dofs.size());
                    for (size_t k = 0; k < emag_dofs.size(); ++k) {
                        nodal_voltages(k) = (emag_dofs[k] != -1) ? emag_solution(emag_dofs[k]) : 0.0;
                    }

                    Eigen::VectorXd nodal_temperatures(heat_dofs.size());
                     for (size_t k = 0; k < heat_dofs.size(); ++k) {
                        nodal_temperatures(k) = (heat_dofs[k] != -1) ? heat_solution(heat_dofs[k]) : 0.0;
                    }

                    double integrated_joule_heat = 0.0;
                    // 假设两个场的FEValues对象使用相同的积分规则
                    for (size_t q_p = 0; q_p < emag_fe_values->num_quadrature_points(); ++q_p) {
                        emag_fe_values->reinit(q_p);
                        heat_fe_values->reinit(q_p);

                        // 获取当前积分点的值
                        const auto& emag_B = emag_fe_values->get_shape_gradients();
                        const auto& heat_N = heat_fe_values->get_shape_values();
                        const double detJ_x_w = emag_fe_values->get_detJ_times_weight();

                        // 1. 在积分点插值计算温度
                        double temp_at_qp = heat_N.transpose() * nodal_temperatures;

                        // 2. 根据插值得到的温度获取电导率
                        const double sigma_at_qp = material.getProperty("electrical_conductivity", {{"Temperature", temp_at_qp}});

                        // 3. 计算积分点的电场强度
                        Eigen::Vector3d grad_V = emag_B * nodal_voltages;

                        // 4. 计算焦耳热密度并累加到积分中
                        integrated_joule_heat += sigma_at_qp * grad_V.squaredNorm() * detJ_x_w;
                    }
                    total_element_power = integrated_joule_heat;
                }
            }
            else if (is_2d) {
                if (auto* tri_elem = dynamic_cast<TriElement*>(elem_ptr)) {
                    tri_elem->setOrder(emag_field_->getElementOrder());
                    auto emag_fe_values = tri_elem->create_fe_values(emag_field_->getElementOrder());

                    tri_elem->setOrder(heat_field_->getElementOrder());
                    auto heat_fe_values = tri_elem->create_fe_values(heat_field_->getElementOrder());

                    const auto emag_dofs = emag_field_->getElementDofs(tri_elem);
                    const auto heat_dofs = heat_field_->getElementDofs(tri_elem);

                    Eigen::VectorXd nodal_voltages(emag_dofs.size());
                    for (size_t k = 0; k < emag_dofs.size(); ++k) {
                       nodal_voltages(k) = (emag_dofs[k] != -1) ? emag_solution(emag_dofs[k]) : 0.0;
                    }

                    Eigen::VectorXd nodal_temperatures(heat_dofs.size());
                    for (size_t k = 0; k < heat_dofs.size(); ++k) {
                       nodal_temperatures(k) = (heat_dofs[k] != -1) ? heat_solution(heat_dofs[k]) : 0.0;
                    }

                    double integrated_joule_heat = 0.0;
                    for (size_t q_p = 0; q_p < emag_fe_values->num_quadrature_points(); ++q_p) {
                        emag_fe_values->reinit(q_p);
                        heat_fe_values->reinit(q_p);

                        const auto& emag_B = emag_fe_values->get_shape_gradients();
                        const auto& heat_N = heat_fe_values->get_shape_values();
                        const double detJ_x_w = emag_fe_values->get_detJ_times_weight();

                        double temp_at_qp = heat_N.transpose() * nodal_temperatures;
                        const double sigma_at_qp = material.getProperty("electrical_conductivity", {{"Temperature", temp_at_qp}});
                        Eigen::Vector2d grad_V = emag_B * nodal_voltages;
                        integrated_joule_heat += sigma_at_qp * grad_V.squaredNorm() * detJ_x_w;
                    }
                    total_element_power = integrated_joule_heat;
                }
            } else if (is_1d) {
                 if (auto* line_elem = dynamic_cast<LineElement*>(elem_ptr)) {
                    line_elem->setOrder(emag_field_->getElementOrder());
                    auto emag_fe_values = line_elem->create_fe_values(emag_field_->getElementOrder());

                    line_elem->setOrder(heat_field_->getElementOrder());
                    auto heat_fe_values = line_elem->create_fe_values(heat_field_->getElementOrder());

                    const auto emag_dofs = emag_field_->getElementDofs(line_elem);
                    const auto heat_dofs = heat_field_->getElementDofs(line_elem);

                    Eigen::VectorXd nodal_voltages(emag_dofs.size());
                    for (size_t k = 0; k < emag_dofs.size(); ++k) {
                        nodal_voltages(k) = (emag_dofs[k] != -1) ? emag_solution(emag_dofs[k]) : 0.0;
                    }

                    Eigen::VectorXd nodal_temperatures(heat_dofs.size());
                    for (size_t k = 0; k < heat_dofs.size(); ++k) {
                        nodal_temperatures(k) = (heat_dofs[k] != -1) ? heat_solution(heat_dofs[k]) : 0.0;
                    }

                    double integrated_joule_heat = 0.0;
                    for (size_t q_p = 0; q_p < emag_fe_values->num_quadrature_points(); ++q_p) {
                        emag_fe_values->reinit(q_p);
                        heat_fe_values->reinit(q_p);

                        const auto& emag_B = emag_fe_values->get_shape_gradients();
                        const auto& heat_N = heat_fe_values->get_shape_values();
                        const double detJ_x_w = emag_fe_values->get_detJ_times_weight();

                        double temp_at_qp = heat_N.transpose() * nodal_temperatures;
                        const double sigma_at_qp = material.getProperty("electrical_conductivity", {{"Temperature", temp_at_qp}});
                        Eigen::Vector<double, 1> grad_V = emag_B * nodal_voltages;
                        integrated_joule_heat += sigma_at_qp * grad_V.squaredNorm() * detJ_x_w;
                    }
                    total_element_power = integrated_joule_heat;
                }
            }

            if (total_element_power > 0) {
                auto source = std::make_unique<VolumetricSource>(
                    elem_ptr->getId(), total_element_power, joule_heat_tag
                );
                heat_field_->addSource(std::move(source));
            }
        }
    }

} // namespace Core