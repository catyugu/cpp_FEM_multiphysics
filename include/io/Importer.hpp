#ifndef IMPORTER_HPP
#define IMPORTER_HPP

#include <map>
#include <string>
#include <memory>
#include <vector>

// Forward declaration
namespace Core {
    class Mesh;
}

namespace IO {

    struct VtuData {
        std::vector<std::vector<double>> points;
        std::map<std::string, std::vector<double>> point_data;
    };

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
        /**
         * @brief Reads nodal data from a VTK Unstructured Grid (.vtu) file.
         * @param filename The path to the input .vtu file.
         * @param data_array_name The name of the DataArray to read (e.g., "Temperature_@_t=0").
         * @return A vector of doubles containing the nodal data.
         */
        static std::vector<double> read_vtu_data(const std::string& filename, const std::string& data_array_name);
        /**
            * @brief Reads node coordinates and specified data arrays from a VTK Unstructured Grid (.vtu) file.
            * @param filename The path to the input .vtu file.
            * @param data_array_names A vector of names of the DataArrays to read.
            * @return A VtuData struct containing the points and a map of the requested data arrays.
            */
        static VtuData read_vtu_points_and_data(const std::string& filename, const std::vector<std::string>& data_array_names);

    };



} // namespace IO

#endif // IMPORTER_HPP
