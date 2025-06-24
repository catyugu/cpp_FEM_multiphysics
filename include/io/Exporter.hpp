#ifndef EXPORTER_HPP
#define EXPORTER_HPP

#include <string>

// Forward declaration
namespace Core {
        class Problem;
}

namespace IO {

        /**
         * @class Exporter
         * @brief A utility class for exporting simulation results to various file formats.
         */
        class Exporter {
        public:
                /**
                 * @brief Writes the mesh and all solved nodal data from a Problem to a legacy VTK file.
                 * @param filename The path to the output file (e.g., "results.vtk").
                 * @param problem The problem instance containing the mesh and solution data.
                 * @return True if the file was written successfully, false otherwise.
                 */
                static bool write_vtk (const std::string& filename, const Core::Problem& problem);
        };

} // namespace IO

#endif // EXPORTER_HPP
