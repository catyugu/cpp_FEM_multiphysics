#include "utils/Profiler.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>

namespace Utils {

    thread_local std::stack<Profiler::ProfileRecord> Profiler::recordStack_;

    Profiler::Profiler() = default;

    Profiler& Profiler::instance() {
        static Profiler instance;
        return instance;
    }

    void Profiler::begin(const std::string& name) {
        if (!enabled_) return;

        ProfileRecord record;
        record.name = name;
        record.startTime = std::chrono::high_resolution_clock::now();
        recordStack_.push(record);
    }

    void Profiler::end() {
        if (!enabled_ || recordStack_.empty()) return;

        auto endTime = std::chrono::high_resolution_clock::now();
        auto record = recordStack_.top();
        recordStack_.pop();

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - record.startTime).count();
        double milliseconds = duration / 1e6;
        std::lock_guard<std::mutex> lock(this->dataMutex_);
        auto& data = this->profileData_[record.name];
        if (!data) {
            data = std::make_unique<ProfileData>();
        }
        data->update(milliseconds);
    }

    std::string Profiler::getReport() const {
        if (!enabled_) return "Profiler is disabled.";

        std::lock_guard<std::mutex> lock(this->dataMutex_);

        if (this->profileData_.empty()) {
            return "No profiling data collected.";
        }

        struct ReportEntry {
            std::string name;
            size_t callCount;
            double totalTime;
            double avgTime;
            double minTime;
            double maxTime;
            double percentage;
        };

        std::vector<ReportEntry> entries;
        double grandTotalTime = 0.0;

        for (const auto& pair : this->profileData_) {
            const auto& name = pair.first;
            const auto& data = pair.second;
            
            std::lock_guard<std::mutex> lock(data->mutex);
            
            if (data->callCount > 0) {
                ReportEntry entry;
                entry.name = name;
                entry.callCount = data->callCount;
                entry.totalTime = data->totalTime;
                entry.avgTime = data->totalTime / data->callCount;
                entry.minTime = data->minTime;
                entry.maxTime = data->maxTime;
                
                entries.push_back(entry);
                grandTotalTime += data->totalTime;
            }
        }

        for (auto& entry : entries) {
            entry.percentage = (entry.totalTime / grandTotalTime) * 100.0;
        }
        std::sort(entries.begin(), entries.end(), [](const ReportEntry& a, const ReportEntry& b) {
            return a.totalTime > b.totalTime;
        });

        std::ostringstream report;
        report << "=== Performance Profiling Report ===\n";
        report << "Total measured time: " << std::fixed << std::setprecision(3) << grandTotalTime << " ms\n";
        report << std::string(60, '-') << "\n";
        report << std::left << std::setw(30) << "Function/Scope" 
               << std::right << std::setw(8) << "Calls" 
               << std::setw(12) << "Total (ms)" 
               << std::setw(9) << "Avg (ms)" 
               << std::setw(9) << "Min (ms)" 
               << std::setw(9) << "Max (ms)" 
               << std::setw(8) << "Pct (%)" << "\n";
        report << std::string(60, '-') << "\n";

        for (const auto& entry : entries) {
            report << std::left << std::setw(30) << entry.name.substr(0, 30)
                   << std::right << std::setw(8) << entry.callCount
                   << std::setw(12) << std::fixed << std::setprecision(3) << entry.totalTime
                   << std::setw(9) << std::fixed << std::setprecision(3) << entry.avgTime
                   << std::setw(9) << std::fixed << std::setprecision(3) << entry.minTime
                   << std::setw(9) << std::fixed << std::setprecision(3) << entry.maxTime
                   << std::setw(7) << std::fixed << std::setprecision(1) << entry.percentage << "%\n";
        }

        return report.str();
    }

    void Profiler::reset() {
        std::lock_guard<std::mutex> lock(this->dataMutex_);
        this->profileData_.clear();
        while (!recordStack_.empty()) {
            recordStack_.pop();
        }
    }

    void Profiler::setEnabled(bool enabled) {
        enabled_ = enabled;
    }

    bool Profiler::isEnabled() const {
        return enabled_;
    }

    ProfileScope::ProfileScope(const std::string& name) {
        Profiler::instance().begin(name);
    }

    ProfileScope::~ProfileScope() {
        Profiler::instance().end();
    }

} // namespace Utils