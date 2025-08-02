#include "post/HeatFluxCalculator.hpp"
#include "core/Problem.hpp"
#include "physics/PhysicsField.hpp"
#include "core/mesh/Element.hpp"
#include "core/FEValues.hpp"
#include "utils/Exceptions.hpp"
#include "utils/SimpleLogger.hpp"

namespace Post {

    const char* HeatFluxCalculator::getName() const {
        return "HeatFlux";
    }

    PostProcessingResult HeatFluxCalculator::compute_derived_quantities(const Core::Problem& problem) const {
        auto& logger = Utils::Logger::instance();
        logger.info("Post-processing: Calculating Heat Flux...");

        const auto* heat_field = problem.getField("Temperature");
        if (!heat_field) {
            throw Exception::ConfigurationException("HeatFluxCalculator requires a 'Temperature' field.");
        }

        const auto& mesh = problem.getMesh();
        const auto& T_solution = heat_field->getSolution();
        const auto& material = heat_field->getMaterial();
        const double k_therm = material.getProperty("thermal_conductivity");

        PostProcessingResult result;
        result.name = getName();
        result.dimension = heat_field->getDimension();
        result.data.resize(mesh.getElements().size());

        for (size_t i = 0; i < mesh.getElements().size(); ++i) {
            Core::Element* elem = mesh.getElements()[i];
            elem->setOrder(heat_field->getElementOrder());

            const auto& ref_data = Core::ReferenceElementCache::get(elem->getTypeName(), elem->getNodes().size(), heat_field->getElementOrder(), heat_field->getElementOrder());
            Core::FEValues fe_values(elem->getGeometry(), heat_field->getElementOrder(), ref_data);

            const auto dofs = heat_field->getElementDofs(elem);

            Eigen::VectorXd nodal_temperatures(dofs.size());
            for (size_t j = 0; j < dofs.size(); ++j) {
                nodal_temperatures(j) = (dofs[j] != -1) ? T_solution(dofs[j]) : 0.0;
            }

            // Store the flux vector for each quadrature point of the element
            result.data[i].resize(fe_values.num_quadrature_points());
            for (size_t q_p = 0; q_p < fe_values.num_quadrature_points(); ++q_p) {
                fe_values.reinit(q_p);
                const auto& B = fe_values.get_shape_gradients(); // B is ∇N
                Eigen::VectorXd grad_T = B * nodal_temperatures;
                result.data[i][q_p] = -k_therm * grad_T;
            }
        }

        logger.info("Heat Flux calculation complete.");
        return result;
    }

} // namespace Post
