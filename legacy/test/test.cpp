#include <gtest/gtest.h>
#include <iostream>
#include <Eigen/Dense>
#include "utils.h"
double SVDTest() {
    Eigen::MatrixXd m_rand = Eigen::MatrixXd::Random(10, 10);
    auto svd_mat = m_rand.jacobiSvd(Eigen::ComputeFullU | Eigen::ComputeFullV);
    const auto& Sigma = svd_mat.singularValues();
    const auto& U = svd_mat.matrixU();
    const auto& V = svd_mat.matrixV();
    // 构造对角矩阵并扩展列数
    Eigen::MatrixXd S = Sigma.asDiagonal();
    // 计算 L^2 Loss
    double L2Loss = (m_rand - U * S * V.transpose()).squaredNorm();
    return L2Loss;
}
TEST(SVDTest, Test) {
    for (int i = 0; i < 10; ++i) {
        double loss = SVDTest();
        EXPECT_LE(loss, 1e-20);
    }
}
TEST(IntegralTest, Test) {
    double res = numerical_integral_1D(
        [](double x) { return x*x; }, 0, 1);
    EXPECT_NEAR(res, 0.333333333333333333333, 1e-10);
    res = numerical_integral_2D(
        [](Eigen::Vector2d x) { return x(0)*x(1); }, 0, 1, 0, 1);
        EXPECT_NEAR(res, 0.25, 1e-10);
    res = numerical_integral_3D(
        [](Eigen::Vector3d x) { return x(0)*x(0)*x(1)*x(2); }, 0, 1, 0, 1, 0, 1);
        EXPECT_NEAR(res, 0.08333333333333333, 1e-10);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}