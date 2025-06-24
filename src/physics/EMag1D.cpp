#include "physics/EMag1D.hpp"
#include "utils/SimpleLogger.hpp"
#include "core/Element.hpp"

namespace Physics {

    EMag1D::EMag1D(double sigma) : mesh_(nullptr), sigma_(sigma) {}

    const char* EMag1D::getName() const {
        return "Electromagnetics 1D";
    }

    void EMag1D::setup(Core::Mesh& mesh) {
        mesh_ = &mesh;
        SimpleLogger::Logger::instance().info("Setting up ", getName(), " for mesh.");
        SimpleLogger::Logger::instance().info("Electrical Conductivity (sigma): ", sigma_);
    }

    void EMag1D::assemble() {
        auto& logger = SimpleLogger::Logger::instance();
        logger.info("Assembling system for ", getName());

        if (!mesh_) {
            logger.error("Mesh is not set up for ", getName());
            return;
        }

        // In a real solver, you would initialize your global matrix and vector here.

        for (const auto& elem_ptr : mesh_->getElements()) {
            auto* line_elem = dynamic_cast<Core::LineElement*>(elem_ptr);
            if(line_elem){
                double h = line_elem->getLength();
                double ke_val = sigma_ / h;

                // Element "stiffness" matrix for 1D EMag: sigma/h * [[1, -1], [-1, 1]]
                // In a real implementation, you would add this to the global matrix.
                logger.info("  - Element ", line_elem->getId(), ", length ", h, ", ke = ", ke_val);
            }
        }
        logger.info("Assembly for ", getName(), " complete.");
    }

} // namespace Physics
