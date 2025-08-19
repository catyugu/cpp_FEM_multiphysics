#ifndef SIMPLELOGGER_HPP
#define SIMPLELOGGER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <mutex>
#include <chrono>
#include <iomanip>

// For enabling colors on Windows
#ifdef _WIN32
#include <windows.h>
#endif

// 前向声明Profiler类，避免循环包含
namespace Utils {
    class Profiler;
}

namespace Utils {
    // Enum for log levels
    enum LogLevel {
        info,
        warn,
        error
    };
    // The Logger class (Singleton)
    class Logger {
    public:
        // Delete copy constructor and assignment operator
        Logger(const Logger &) = delete;

        Logger &operator=(const Logger &) = delete;

        // Get the single instance of the logger
        static Logger &instance() {
            static Logger instance;
            return instance;
        }

        // Set log file, optionally with a max size
        void set_logfile(const std::string &log_filename = "") {
            if (!log_filename.empty()) {
                this->log_file_name = log_filename;
                log_file.open(log_filename, std::ios::out | std::ios::app);
                if (!log_file.is_open()) {
                    log(LogLevel::error, __FILE__,"Failed to open log file: ", log_filename);
                }
            }
            if (log_file.is_open() && log_filename.empty()) {
                log_file.close();
            }
        }
        void set_loglevel(int log_level = 0) {
            if (log_level < 0 || log_level > 2) {
                log(LogLevel::warn, __FILE__, "Invalid log level: ", log_level);
                log(LogLevel::warn, __FILE__, "Using default log level: ", 0);
                this->log_level = 0;
                return;
            }
            this->log_level = log_level;
        }

        // Variadic template log functions for different levels
        template<typename... Args>
        void info(Args... args) {
            log(LogLevel::info, args...);
        }

        template<typename... Args>
        void warn(Args... args) {
            log(LogLevel::warn, args...);
        }

        template<typename... Args>
        void error(Args... args) {
            log(LogLevel::error, args...);
        }

    private:
        Logger() {
#ifdef _WIN32
            enable_windows_virtual_terminal();
#endif
        }

        // Private destructor to close the file
        ~Logger() {
            if (log_file.is_open()) {
                log_file.close();
            }
        }

        // The core logging function
        template<typename... Args>
        void log(LogLevel level,
            Args... args) {
            if (level < log_level) {
                return;
            }
            // Use a stringstream to build the message from variadic arguments
            std::stringstream message_stream;
            // This is a C++17 fold expression to stream all arguments
            ((message_stream << args), ...);

            std::string message = message_stream.str();

            // Lock the mutex to ensure thread safety
            std::lock_guard<std::mutex> lock(log_mutex);
            // Get current time
            auto now = std::chrono::system_clock::now();
            auto in_time_t = std::chrono::system_clock::to_time_t(now);

            std::stringstream timestamp_stream;
#if defined(_MSC_VER)
            struct tm timeinfo{};
            localtime_s(&timeinfo, &in_time_t);
            timestamp_stream << std::put_time(&timeinfo, "%Y-%m-%d %X");
#else
            timestamp_stream << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
#endif
            // Log to console with colors
            std::cout << get_level_color(level)
                    << "[" << get_level_string(level) << "] "
                    << "\033[0m"
                    << "[" << timestamp_stream.str() << "] "
                    << message
                    << std::endl;

            // Log to file if it's open
            if (log_file.is_open()) {
                log_file << "[" << timestamp_stream.str() << "] "
                        << "[" << get_level_string(level) << "] "
                        << message
                        << std::endl;
            }

        }

        // Helper to get the string representation of a level
        static std::string get_level_string(LogLevel level) {
            switch (level) {
                case LogLevel::info: return "INFO";
                case LogLevel::warn: return "WARN";
                case LogLevel::error: return "ERROR";
                default: return "?????";
            }
        }

        // Helper to get the ANSI color code for a level
        static const char *get_level_color(LogLevel level) {
            switch (level) {
                case LogLevel::info: return "\033[32m"; // Green
                case LogLevel::warn: return "\033[1;33m"; // Bold Yellow
                case LogLevel::error: return "\033[1;31m"; // Bold Red
                default: return "\033[0m"; // Reset
            }
        }

#ifdef _WIN32
        // Helper to enable virtual terminal processing on Windows
        static void enable_windows_virtual_terminal() {
            HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
            if (hOut == INVALID_HANDLE_VALUE) return;
            DWORD dwMode = 0;
            if (!GetConsoleMode(hOut, &dwMode)) return;
            dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
            SetConsoleMode(hOut, dwMode);
        }
#endif

        std::mutex log_mutex;
        std::ofstream log_file;
        std::string log_file_name = "";
        int log_level = 0;
    };
} // namespace Utils

#endif // SIMPLELOGGER_HPP