// catyugu/cpp_fem_multiphysics/cpp_FEM_multiphysics-dev/tests/test_fe_values.cpp
#include <gtest/gtest.h>
#include <utils/SimpleLogger.hpp>

#include "core/mesh/TriElement.hpp"
#include "core/mesh/Node.hpp"
#include "core/FEValues.hpp"

TEST(FEValuesTest, TriangleAreaAndGradients) {
    auto& logger = Utils::Logger::instance();
    // Create a simple right-angled triangle
    auto n1 = new Core::Node(0, 0, 0);
    auto n2 = new Core::Node(1, 1, 0);
    auto n3 = new Core::Node(2, 0, 1);
    auto tri = new Core::TriElement(0);
    tri->addNode(n1);
    tri->addNode(n2);
    tri->addNode(n3);
    tri->setOrder(1); // Linear element
    logger.info("FEValuesTest: TriangleAreaAndGradients");
    logger.info("FEValuesTest: TriangleAreaAndGradients: Area: ", tri->getArea());

    // Create FEValues with a 1-point quadrature rule (exact for linear triangles)
    auto fe_values = tri->create_fe_values(1);

    logger.info("FEValuesTest: TriangleAreaAndGradients: Quadrature points: ", fe_values->num_quadrature_points());



    // There should be one quadrature point
    ASSERT_EQ(fe_values->num_quadrature_points(), 1);

    // The area of this triangle is 0.5.
    // For a 1-point rule, the weight is 0.5, so det(J) should be 1.0.
    // And det(J) * weight should be 0.5.
    ASSERT_NEAR(fe_values->get_detJ_times_weight(), 0.5, 1e-9);

    // For a linear triangle, the gradients (B-matrix) should be constant.
    const auto& B = fe_values->get_shape_gradients();
    // B = (1/2A) * [[y23, y31, y12], [x32, x13, x21]]
    // A = 0.5, so B = [[-1, 1, 0], [-1, 0, 1]]
    Eigen::MatrixXd B_expected(2, 3);
    B_expected << -1.0, 1.0, 0.0,
                  -1.0, 0.0, 1.0;

    ASSERT_TRUE(B.isApprox(B_expected, 1e-9));

    delete tri;
    delete n1;
    delete n2;
    delete n3;
}