#ifndef IMPORTER_HPP
#define IMPORTER_HPP

#include <string>
#include <memory>

// Forward declaration
namespace Core {
    class Mesh;
}

namespace IO {

    /**
     * @class Importer
     * @brief A utility class for importing mesh data from various file formats.
     */
    class Importer {
    public:
        /**
         * @brief Reads a mesh from a COMSOL .mphtxt text file.
         * This parser expects a specific format containing vertex coordinates and triangular elements.
         * @param filename The path to the input .mphtxt file.
         * @return A unique_ptr to the newly created Mesh object, or nullptr on failure.
         */
        static std::unique_ptr<Core::Mesh> read_comsol_mphtxt(const std::string& filename);
        /**
         * @brief Reads a mesh from a GMSH .msh file.
         * This parser expects a specific format containing vertex coordinates and triangular elements.
         * @param filename The path to the input .msh file.
         * @return A unique_ptr to the newly created Mesh object, or nullptr on failure.
         */
        std::unique_ptr<Core::Mesh> Importer::read_gmsh_msh(const std::string& filename);
    };

} // namespace IO

#endif // IMPORTER_HPP
