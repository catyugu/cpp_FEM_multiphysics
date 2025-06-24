#include "physics/Heat1D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Element.hpp"


namespace Physics {

    Heat1D::Heat1D(double k) : mesh_(nullptr), k_(k) {}

    const char* Heat1D::getName() const {
        return "Heat Transfer 1D";
    }

    void Heat1D::setup(Core::Mesh& mesh) {
        mesh_ = &mesh;
        SimpleLogger::Logger::instance().info("Setting up ", getName(), " for mesh.");
        SimpleLogger::Logger::instance().info("Thermal Conductivity (k): ", k_);
    }

    void Heat1D::assemble() {
        auto& logger = SimpleLogger::Logger::instance();
        logger.info("Assembling system for ", getName());

        if (!mesh_) {
            logger.error("Mesh is not set up for ", getName());
            return;
        }

        for (const auto& elem_ptr : mesh_->getElements()) {
            auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
            if(line_elem){
                double h = line_elem->getLength();
                double ke_val = k_ / h;

                // Element stiffness matrix for 1D Heat: k/h * [[1, -1], [-1, 1]]
                logger.info("  - Element ", line_elem->getId(), ", length ", h, ", ke = ", ke_val);

                // Element force vector from heat source Q would be calculated here
                // Q * h/2 * {1, 1}
                // This is where the coupling happens. The 'heat_source_' vector
                // would be used to calculate the element force vector.
            }
        }

        logger.info("Assembly for ", getName(), " complete.");
    }

    void Heat1D::setHeatSource(const std::vector<double>& source) {
        heat_source_ = source;
        SimpleLogger::Logger::instance().info(getName(), ": Received heat source data.");
    }


} // namespace Physics
