#ifndef REFERENCEELEMENT_HPP
#define REFERENCEELEMENT_HPP

#include "utils/Quadrature.hpp"
#include "utils/ShapeFunctions.hpp"
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <tuple>
#include <map>
#include <mutex>

namespace Core {

    // 存储预计算的参考单元数据
    struct ReferenceElementData {
        std::vector<Utils::QuadraturePoint> quadrature_points;
        std::vector<Eigen::VectorXd> N_values_at_qps;         // 所有积分点上的形函数值
        std::vector<Eigen::MatrixXd> dN_d_natural_at_qps; // 所有积分点上，在自然坐标系下的导数
    };

    // 用于管理参考单元数据缓存的单例类
    class ReferenceElementCache {
    public:
        // 使用元组作为缓存的键: {单元类型名, 单元阶次, 积分阶次}
        using CacheKey = std::tuple<std::string, int, int>;

        /**
         * @brief 获取指定类型的参考单元数据。如果缓存中不存在，则会先计算再存入缓存。
         * @param type_name 单元类型名 (e.g., "TriElement")
         * @param num_vertices 单元的几何顶点数
         * @param fe_order 有限元近似的阶次
         * @param quad_order 积分规则的阶次
         * @return 对缓存中参考单元数据的常量引用
         */
        static const ReferenceElementData& get(const std::string& type_name, int num_vertices, int fe_order, int quad_order);

    private:
        static std::map<CacheKey, ReferenceElementData> cache_;
        static std::mutex mutex_;
    };

} // namespace Core
#endif // REFERENCEELEMENT_HPP