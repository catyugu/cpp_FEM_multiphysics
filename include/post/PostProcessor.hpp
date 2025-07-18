#ifndef POSTPROCESSOR_HPP
#define POSTPROCESSOR_HPP

#include <string>
#include <vector>
#include <Eigen/Dense>

// Forward declaration
namespace Core {
    class Problem;
}

namespace Post {

    /**
     * @struct PostProcessingResult
     * @brief Holds the results of a post-processing computation.
     *
     * This struct stores derived quantities, such as heat flux or stress,
     * which are typically calculated per element.
     */
    struct PostProcessingResult {
        std::string name;
        int dimension; // 1 for scalar, 3 for vector, etc.
        // The outer vector corresponds to the element index.
        // The inner vector holds the computed values for that element (e.g., at quadrature points).
        std::vector<std::vector<Eigen::VectorXd>> data;
    };

    /**
     * @class PostProcessor
     * @brief Abstract base class for all post-processing computations.
     *
     * This class defines the interface for calculating derived quantities
     * from the primary solution variables after a simulation has run.
     */
    class PostProcessor {
    public:
        virtual ~PostProcessor() = default;

        /**
         * @brief Computes derived quantities from the problem's solution.
         * @param problem The problem instance containing the mesh and solution data.
         * @return A PostProcessingResult struct containing the computed data.
         */
        virtual PostProcessingResult compute_derived_quantities(const Core::Problem& problem) const = 0;

        /**
         * @brief Gets the name of the derived quantity this post-processor calculates.
         * @return The name of the result (e.g., "HeatFlux").
         */
        virtual const char* getName() const = 0;
    };

} // namespace Post

#endif // POSTPROCESSOR_HPP
