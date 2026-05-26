//
// Created by Diaz, Diego on 9.12.2024.
//
#ifndef CDS_LOGGER_H
#define CDS_LOGGER_H

#include <iostream>
#include <string>
#include "utils.h"
#ifdef USE_MALLOC_COUNT
#include <malloc_count.h>
#endif

//code generated with ChatGPT
#define CONCAT_IMPL(x,y) x##y
#define CONCAT(x,y) CONCAT_IMPL(x,y)

// -----------------------------
// Log levels
// -----------------------------
enum class log_level {
    ERROR = 0,
    WARN  = 1,
    INFO  = 2,
    DEBUG = 3,
    TRACE = 4
};

// -----------------------------
// Logger
// -----------------------------
struct Logger {

    static log_level level;
    static std::chrono::steady_clock::time_point start;
    static int indent;
    //static std::chrono::steady_clock::time_point last;

    static bool enabled(log_level msgLevel) {
        return static_cast<int>(msgLevel) <= static_cast<int>(level);
    }

    static void log(log_level msgLevel, const std::string& msg) {
        if (!enabled(msgLevel))
            return;

        auto now = std::chrono::steady_clock::now();

        double total = std::chrono::duration<double>(now - start).count();
        //double delta = std::chrono::duration<double>(now - last).count();
        //last = now;

        int hours   = static_cast<int>(total) / 3600;
        int minutes = (static_cast<int>(total) % 3600) / 60;
        int seconds = static_cast<int>(total) % 60;

        std::cout << "["
                  << std::setw(2) << std::setfill('0') << hours << ":"
                  << std::setw(2) << minutes << ":"
                  << std::setw(2) << seconds << "] ";

        for (int i = 0; i < indent; ++i)
            std::cout << "  ";
        std::cout<<msg<<std::endl;
        // Only show delta for TRACE
        //if (msgLevel >= log_level::TRACE) std::cout << " | +" << format_delta(delta);
        //std::cout << "] " << msg << std::endl;
    }
};

int Logger::indent = 0;
log_level Logger::level = log_level::INFO;
std::chrono::steady_clock::time_point Logger::start = std::chrono::steady_clock::now();
//std::chrono::steady_clock::time_point Logger::last = start;

class ScopeIndent {
public:
    ScopeIndent() {
        Logger::indent++;
    }
    ~ScopeIndent() {
        Logger::indent--;
    }
};

// -----------------------------
// RAII Scope Logger
// -----------------------------
class ScopeLog {
public:
    ScopeLog(const std::string& name, log_level lvl)
        : name(name),
          level(lvl),
          active(Logger::enabled(lvl)),
          start(std::chrono::steady_clock::now()) {
        if (active) {
            Logger::log(level, "→ " + name);
#if USE_MALLOC_COUNT
            malloc_count_reset_peak();
#endif
        }
    }

    ~ScopeLog() {
        if (!active) return;
        auto end = std::chrono::steady_clock::now();

        double elapsed = std::chrono::duration<double>(end - start).count();
        std::string msg = "← " + name + " (" + format_time(elapsed);
#if USE_MALLOC_COUNT
        msg += ", " + format_space(static_cast<off_t>(malloc_count_peak()));
#endif
        msg += ")";
        Logger::log(level, msg);
    }

private:
    std::string name;
    log_level level;
    bool active;
    std::chrono::steady_clock::time_point start;
};

// -----------------------------
// Logging macros
// -----------------------------
#define TRACE_SCOPE() ScopeLog scope(__FUNCTION__, log_level::TRACE)
#define DEBUG_SCOPE() ScopeLog scope(__FUNCTION__, log_level::DEBUG)
#define SCOPE_INFO() ScopeIndent CONCAT(_scopeIndent_, __LINE__)

#define LOG_ERROR(msg) Logger::log(log_level::ERROR, msg)
#define LOG_WARN(msg)  Logger::log(log_level::WARN,  msg)
#define LOG_INFO(msg)  Logger::log(log_level::INFO,  msg)
#define LOG_DEBUG(msg) Logger::log(log_level::DEBUG, msg)
#define LOG_TRACE(msg) Logger::log(log_level::TRACE, msg)
#endif //CDS_LOGGER_H
