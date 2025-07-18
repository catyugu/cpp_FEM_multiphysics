#ifndef HEATFLUXCALCULATOR_HPP
#define HEATFLUXCALCULATOR_HPP

#include "PostProcessor.hpp"

namespace Post {

    /**
     * @class HeatFluxCalculator
     * @brief A concrete post-processor for calculating heat flux from a temperature solution.
     *
     * This class implements the logic to compute the heat flux vector (q = -k∇T)
     * for each element in the mesh based on the solved temperature field.
     */
    class HeatFluxCalculator : public PostProcessor {
    public:
        PostProcessingResult compute_derived_quantities(const Core::Problem& problem) const override;
        const char* getName() const override;
    };

} // namespace Post

#endif // HEATFLUXCALCULATOR_HPP
