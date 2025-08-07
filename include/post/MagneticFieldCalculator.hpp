#ifndef MAGNETICFIELDCALCULATOR_HPP
#define MAGNETICFIELDCALCULATOR_HPP

#include "PostProcessor.hpp"

namespace Post {
    class MagneticFieldCalculator : public PostProcessor {
    public:
        PostProcessingResult compute_derived_quantities(const Core::Problem &problem) const override;

        const char *getName() const override;
    };
} // namespace Post

#endif // MAGNETICFIELDCALCULATOR_HPP
