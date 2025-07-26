#include "core/ReferenceElement.hpp"
#include "utils/Exceptions.hpp"
#include <mutex>

namespace Core {

// 初始化静态成员
std::map<ReferenceElementCache::CacheKey, ReferenceElementData> ReferenceElementCache::cache_;
std::mutex ReferenceElementCache::mutex_;

const ReferenceElementData& ReferenceElementCache::get(const std::string& type_name, int num_vertices, int fe_order, int quad_order) {
    CacheKey key = {type_name, fe_order, quad_order};

    // 使用互斥锁保证线程安全
    std::lock_guard<std::mutex> lock(mutex_);
    
    // 如果缓存中未找到，则计算并存入
    if (cache_.find(key) == cache_.end()) {
        ReferenceElementData data;

        // 这部分逻辑是从旧的 FEValues 构造函数中抽离出来的、与几何无关的计算
        if (type_name == "LineElement") {
            data.quadrature_points = Utils::Quadrature::getLineQuadrature(quad_order);
            for (const auto& qp : data.quadrature_points) {
                data.N_values_at_qps.push_back(Utils::ShapeFunctions::getLineShapeFunctions(fe_order, qp.point(0)));
                data.dN_d_natural_at_qps.push_back(Utils::ShapeFunctions::getLineShapeFunctionDerivatives(fe_order, qp.point(0)));
            }
        } else if (type_name == "TriElement") {
            data.quadrature_points = Utils::Quadrature::getTriangleQuadrature(quad_order);
            for (const auto& qp : data.quadrature_points) {
                data.N_values_at_qps.push_back(Utils::ShapeFunctions::getTriShapeFunctions(fe_order, qp.point(0), qp.point(1)));
                data.dN_d_natural_at_qps.push_back(Utils::ShapeFunctions::getTriShapeFunctionDerivatives(fe_order, qp.point(0), qp.point(1)));
            }
        } else if (type_name == "TetElement") {
            data.quadrature_points = Utils::Quadrature::getTetrahedronQuadrature(quad_order);
            for (const auto& qp : data.quadrature_points) {
                data.N_values_at_qps.push_back(Utils::ShapeFunctions::getTetShapeFunctions(fe_order, qp.point(0), qp.point(1), qp.point(2)));
                data.dN_d_natural_at_qps.push_back(Utils::ShapeFunctions::getTetShapeFunctionDerivatives(fe_order, qp.point(0), qp.point(1), qp.point(2)));
            }
        } else {
            throw Exception::ConfigurationException("Unsupported element type for ReferenceElementCache: " + type_name);
        }

        cache_[key] = std::move(data);
    }
    
    return cache_.at(key);
}

} // namespace Core