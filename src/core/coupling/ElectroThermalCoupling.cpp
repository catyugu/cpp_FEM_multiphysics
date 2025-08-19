#include "core/coupling/ElectroThermalCoupling.hpp"
#include "core/mesh/Element.hpp"
#include "physics/PhysicsField.hpp"
#include "utils/SimpleLogger.hpp"
#include <string>
#include <core/FEValues.hpp>
#include "core/ReferenceElement.hpp"
#include <omp.h> // 引入OpenMP头文件

namespace Core {
    // ... setup函数保持不变 ...
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
        const auto& elements = mesh->getElements();
        if (elements.empty()) return;

        const auto& emag_solution = emag_field_->getSolution();
        const auto& heat_solution = heat_field_->getSolution();

        Eigen::VectorXd& coupling_rhs = heat_field_->getCouplingRHS();
        coupling_rhs.setZero();

        Eigen::VectorXd global_coupling_rhs = Eigen::VectorXd::Zero(coupling_rhs.size());

        #pragma omp parallel
        {
            Eigen::VectorXd local_coupling_rhs = Eigen::VectorXd::Zero(coupling_rhs.size());

            // --- 修正点：将范围for循环改为传统的索引for循环 ---
            #pragma omp for nowait schedule(static)
            for (int i = 0; i < elements.size(); ++i) {
                const auto& elem_ptr = elements[i]; // 在循环内部获取元素指针

                const auto emag_dofs = emag_field_->getElementDofs(elem_ptr);
                const auto heat_dofs = heat_field_->getElementDofs(elem_ptr);

                Eigen::VectorXd nodal_voltages(static_cast<Eigen::Index>(emag_dofs.size()));
                for (size_t j = 0; j < emag_dofs.size(); ++j) {
                    nodal_voltages(static_cast<Eigen::Index>(j)) = (emag_dofs[j] != -1) ? emag_solution(emag_dofs[j]) : 0.0;
                }

                Eigen::VectorXd nodal_temperatures(static_cast<Eigen::Index>(heat_dofs.size()));
                for (size_t j = 0; j < heat_dofs.size(); ++j) {
                    nodal_temperatures(static_cast<Eigen::Index>(j)) = (heat_dofs[j] != -1) ? heat_solution(heat_dofs[j]) : 0.0;
                }

                // 使用Element内部缓存的FEValues对象，提高性能
                FEValues* emag_fe_values = elem_ptr->getFEValues(2 * emag_field_->getElementOrder(), Core::AnalysisType::SCALAR_DIFFUSION);
                FEValues* heat_fe_values = elem_ptr->getFEValues(2 * heat_field_->getElementOrder(), Core::AnalysisType::SCALAR_DIFFUSION);

                const auto& material = emag_field_->getMaterial(elem_ptr);

                Eigen::VectorXd fe_local = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(heat_dofs.size()));

                for (size_t q_p = 0; q_p < emag_fe_values->num_quadrature_points(); ++q_p) {
                    emag_fe_values->reinit(static_cast<int>(q_p));
                    heat_fe_values->reinit(static_cast<int>(q_p));

                    const auto& emag_B = emag_fe_values->get_shape_gradients();
                    const auto& heat_N = heat_fe_values->get_shape_values();
                    const double detJ_x_w = emag_fe_values->get_detJ_times_weight();

                    double temp_at_qp = heat_N.transpose() * nodal_temperatures;
                    const double sigma_at_qp = material.getProperty("electrical_conductivity", {{"Temperature", temp_at_qp}});

                    Eigen::VectorXd grad_V = emag_B * nodal_voltages;
                    double joule_heating = sigma_at_qp * grad_V.squaredNorm();

                    fe_local += heat_N * joule_heating * detJ_x_w;
                }

                for (size_t j = 0; j < heat_dofs.size(); ++j) {
                    if (heat_dofs[j] != -1) {
                        local_coupling_rhs(heat_dofs[j]) += fe_local(j);
                    }
                }
            }

            #pragma omp critical
            {
                global_coupling_rhs += local_coupling_rhs;
            }
        }

        coupling_rhs = global_coupling_rhs;
    }

} // namespace Core