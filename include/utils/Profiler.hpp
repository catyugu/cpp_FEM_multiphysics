#ifndef PROFILER_HPP
#define PROFILER_HPP

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <stack>
#include <thread>
#include <mutex>
#include <limits>
#include <algorithm>

namespace Utils {
    /**
     * @class Profiler
     * @brief 高度智能化的性能分析器，用于分析代码性能瓶颈
     */
    class Profiler {
    public:
        static Profiler &instance();

        Profiler(const Profiler &) = delete;

        Profiler &operator=(const Profiler &) = delete;

        /**
         * @brief 开始测量一个函数或代码段的性能
         * @param name 测量点的名称
         */
        void begin(const std::string &name);

        /**
         * @brief 结束当前测量点的性能测量
         */
        void end();

        /**
         * @brief 获取性能报告
         * @return 格式化的性能分析报告字符串
         */
        std::string getReport() const;

        /**
         * @brief 重置所有统计数据
         */
        void reset();

        /**
         * @brief 启用或禁用性能分析
         * @param enabled true启用，false禁用
         */
        void setEnabled(bool enabled);

        /**
         * @brief 检查性能分析是否启用
         * @return true如果启用，false如果禁用
         */
        bool isEnabled() const;

    private:
        Profiler();

        struct ProfileData {
            size_t callCount = 0;
            double totalTime = 0.0;
            double minTime = (std::numeric_limits<double>::max)();
            double maxTime = 0.0;
            mutable std::mutex mutex;

            void update(double duration) {
                std::lock_guard<std::mutex> lock(mutex);
                callCount++;
                totalTime += duration;
                if (duration < minTime) minTime = duration;
                if (duration > maxTime) maxTime = duration;
            }
        };

        struct ProfileRecord {
            std::string name;
            std::chrono::high_resolution_clock::time_point startTime;
        };

        static thread_local std::stack<ProfileRecord> recordStack_;

        mutable std::unordered_map<std::string, std::unique_ptr<ProfileData> > profileData_;
        mutable std::mutex dataMutex_;

        bool enabled_ = true;
    };

    /**
     * @class ProfileScope
     * @brief RAII风格的性能分析范围类
     * 
     * 在作用域开始时自动开始性能测量，在作用域结束时自动结束测量。
     */
    class ProfileScope {
    public:
        /**
         * @brief 构造函数，自动开始性能测量
         * @param name 测量点名称
         */
        explicit ProfileScope(const std::string &name);

        /**
         * @brief 析构函数，自动结束性能测量
         */
        ~ProfileScope();

    private:
        ProfileScope(const ProfileScope &) = delete;

        ProfileScope &operator=(const ProfileScope &) = delete;
    };
} // namespace Utils

// 方便使用的宏定义
#define PROFILE_FUNCTION() Utils::ProfileScope profile_scope_##__LINE__(__FUNCTION__)
#define PROFILE_SCOPE(name) Utils::ProfileScope profile_scope_##__LINE__(name)

#endif // PROFILER_HPP
