#include <core/mesh/Element.hpp>
#include <core/FEValues.hpp>
#include <stdexcept>
#include <algorithm>

namespace Core {

    Element::Element(int id) : id_(id), order_(1), geometry_(nullptr) {}

    int Element::getId() const {
        return id_;
    }

    const std::vector<Node*>& Element::getNodes() const {
        return nodes_;
    }

    void Element::addNode(Node* node) {
        nodes_.push_back(node);
    }

    void Element::set_nodes_internal(const std::vector<Node*>& new_nodes) {
        nodes_.clear();
        nodes_.insert(nodes_.begin(), new_nodes.begin(), new_nodes.end());
    }

    void Element::update_geometry() {
        if (!nodes_.empty()) {
            geometry_ = std::make_unique<ElementGeometry>(nodes_, this->getDimension());
        }
    }

    // 新增 getGeometry 的实现
    const ElementGeometry& Element::getGeometry() {
        if (!geometry_) {
            update_geometry();
        }
        return *geometry_;
    }

    // 新增：获取与此单元关联的 FEValues 对象，如果尚未创建则创建并初始化
    FEValues* Element::getFEValues(int quad_order, AnalysisType analysis_type) {
        // 检查缓存是否需要更新 (例如，如果积分阶数或分析类型改变)
        if (!fe_values_ || stored_quad_order_ != quad_order || stored_analysis_type_ != analysis_type) {
} // namespace Core