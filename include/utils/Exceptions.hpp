#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <stdexcept>
#include <string>

namespace Core {

    class FileIOException : public std::runtime_error {
    public:
        explicit FileIOException(const std::string& message)
            : std::runtime_error("File I/O Error: " + message) {}
    };

    class SolverException : public std::runtime_error {
    public:
        explicit SolverException(const std::string& message)
            : std::runtime_error("Solver Error: " + message) {}
    };

    class ConfigurationException : public std::runtime_error {
    public:
        explicit ConfigurationException(const std::string& message)
            : std::runtime_error("Configuration Error: " + message) {}
    };

} // namespace Core

#endif // EXCEPTIONS_HPP