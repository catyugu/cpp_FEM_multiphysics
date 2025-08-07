#ifndef ANALYSIS_TYPES_HPP
#define ANALYSIS_TYPES_HPP

namespace Core {

/**
 * @enum AnalysisType
 * @brief 定义不同的物理分析类型，用于确定B矩阵的构建方式
 */
enum class AnalysisType {
    SCALAR_DIFFUSION,    // 标量扩散问题 (热传导、电流等)，B = ∇N
    VECTOR_CURL,         // 矢量旋度问题 (磁场)，B = curl(N)
    VECTOR_GRADIENT,     // 矢量梯度问题 (固体力学)，B = 应变位移矩阵
    CUSTOM              // 自定义类型，需要特殊处理
};

} // namespace Core

#endif // ANALYSIS_TYPES_HPP
