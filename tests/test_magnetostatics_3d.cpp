#include <gtest/gtest.h>
#include <memory>
#include "core/Problem.hpp"
#include "core/Material.hpp"
#include "physics/Magnetic3D.hpp"
#include "core/sources/PrescribedCurrentDensity.hpp"
#include <core/bcs/BoundaryCondition.hpp>
#include "utils/SimpleLogger.hpp"
#include "post/MagneticFieldCalculator.hpp"

TEST(Magnetostatics3DTest, SolenoidQuantitativeValidation) {
    auto& logger = Utils::Logger::instance();
    logger.info("--- Setting up 3D Magnetostatics Test: Solenoid ---");

    auto mesh = std::unique_ptr<Core::Mesh>(Core::Mesh::create_uniform_3d_mesh(2.0, 2.0, 10.0, 10, 10, 20));
    ASSERT_NE(mesh, nullptr);

    Core::Material air("Air");
    const double mu_0 = 1.25663706e-6;
    air.setProperty("magnetic_permeability", mu_0);

    auto problem = std::make_unique<Core::Problem>(std::move(mesh));
    auto magnetic_field_ptr = std::make_unique<Physics::Magnetic3D>(air);
    problem->addField(std::move(magnetic_field_ptr));
    problem->addPostProcessor(std::make_unique<Post::MagneticFieldCalculator>());
    problem->setup();

    auto* field = problem->getField("MagneticVectorPotential");
    ASSERT_NE(field, nullptr);

    const double current_density_magnitude = 1.0e6;

    for (const auto& elem : problem->getMesh().getElements()) {
        std::vector<double> centroid = {0.0, 0.0, 0.0};
        const auto& elem_nodes = elem->getNodes();
        for (const auto* node : elem_nodes) {
            const auto& coords = node->getCoords();
            centroid[0] += coords[0];
            centroid[1] += coords[1];
        }
        centroid[0] /= elem_nodes.size();
        centroid[1] /= elem_nodes.size();

        double x = centroid[0] - 1.0; // 将坐标中心移至 (0,0)
        double y = centroid[1] - 1.0;
        double r = std::sqrt(x*x + y*y);

        if (r > 0.2 && r < 0.3) {
            Eigen::Vector3d J(-y/r * current_density_magnitude, x/r * current_density_magnitude, 0);
            field->addSource(std::make_unique<Core::PrescribedCurrentDensity>(elem->getId(), J));
        }
    }
    auto& dof_manager = problem->getDofManager();
    for (const auto& node : problem->getMesh().getNodes()) {
        const auto& coords = node->getCoords();
        if (std::abs(coords[0] - 0.0) < 1e-9 || std::abs(coords[0] - 2.0) < 1e-9 ||
            std::abs(coords[1] - 0.0) < 1e-9 || std::abs(coords[1] - 2.0) < 1e-9 ||
            std::abs(coords[2] - 0.0) < 1e-9 || std::abs(coords[2] - 10.0) < 1e-9)
        {
            int base_dof_idx = dof_manager.getEquationIndex(node->getId(), "MagneticVectorPotential");
            if (base_dof_idx != -1) {
                field->addBC(std::make_unique<Core::DirichletBC>(base_dof_idx + 0, Eigen::Vector<double, 1>(0.0)));
                field->addBC(std::make_unique<Core::DirichletBC>(base_dof_idx + 1, Eigen::Vector<double, 1>(0.0)));
                field->addBC(std::make_unique<Core::DirichletBC>(base_dof_idx + 2, Eigen::Vector<double, 1>(0.0)));
            }
        }
    }

    // 求解之前，我们需要先在 `Magnetic3D::assemble` 中修正一个逻辑错误
    ASSERT_NO_THROW(problem->solveSteadyState());
    problem->exportResults("solenoid_results_validated.vtk");

    const auto& post_results = problem->getPostProcessingResults();
    ASSERT_TRUE(post_results.count("MagneticField"));
    const auto& b_field_results = post_results.at("MagneticField").data;

    const double coil_thickness = 0.1;
    const double B_z_analytical = mu_0 * current_density_magnitude * coil_thickness;
    logger.info("Analytical B_z field should be approximately: ", B_z_analytical);

    int points_validated = 0;
    for (size_t i = 0; i < problem->getMesh().getElements().size(); ++i) {
        auto* elem = problem->getMesh().getElement(i);

        // **【修正】**：计算并使用单元的质心进行判断
        std::vector<double> centroid = {0.0, 0.0, 0.0};
        const auto& elem_nodes = elem->getNodes();
        for (const auto* node : elem_nodes) {
            const auto& coords = node->getCoords();
            centroid[0] += coords[0];
            centroid[1] += coords[1];
            centroid[2] += coords[2];
        }
        centroid[0] /= elem_nodes.size();
        centroid[1] /= elem_nodes.size();
        centroid[2] /= elem_nodes.size();

        double x = centroid[0] - 1.0;
        double y = centroid[1] - 1.0;
        double z = centroid[2];
        double r = std::sqrt(x*x + y*y);

        if (r < 0.1 && z > 2.0 && z < 8.0) {
            const auto& B_at_qp0 = b_field_results[i][0];

            ASSERT_NEAR(B_at_qp0(0), 0.0, B_z_analytical * 0.05);
            ASSERT_NEAR(B_at_qp0(1), 0.0, B_z_analytical * 0.05);
            ASSERT_NEAR(B_at_qp0(2), B_z_analytical, B_z_analytical * 0.1);
            points_validated++;
        }
    }
    ASSERT_GT(points_validated, 0);
    logger.info("Validated ", points_validated, " points inside the solenoid against the analytical solution.");
}