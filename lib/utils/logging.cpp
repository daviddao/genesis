/**
 * @brief Implementation of Logging functions.
 *
 * @file
 * @ingroup utils
 */

#include "utils/logging.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#ifdef PTHREADS
#    include <mutex>
#endif

#include "utils/utils.hpp"

namespace genesis {

// =============================================================================
//     Settings
// =============================================================================

#ifdef PTHREADS
    static std::mutex log_mutex;
#endif

// TODO use different init for log details depending on DEBUG

// init static members
LoggingDetails Logging::details = {
    false, // count
    false, // date
    false, // time
    false, // runtime
    false, // rundiff
    false, // file
    false, // line
    false, // function
    true   // level
};
Logging::LoggingLevel      Logging::max_level_  = kDebug4;
long                       Logging::count_      = 0;
clock_t                    Logging::last_clock_ = 0;
std::vector<std::ostream*> Logging::ostreams_;
int                        Logging::report_percentage_ = 5;
std::string                Logging::debug_indent       = "    ";

/**
 * @brief Set the highest log level that is reported.
 *
 * Invocations of log with higher levels will create no output.
 * It creates a warning if the set level is higher than the static compile time
 * level set by #LOG_LEVEL_MAX.
 */
void Logging::max_level (const LoggingLevel level)
{
    if (level > LOG_LEVEL_MAX) {
        LOG_WARN << "Logging max level set to " << level << ", but compile "
                 << "time max level is " << LOG_LEVEL_MAX << ", so that "
                 << "everything above that will not be logged.";
    }
    max_level_ = level;
}

/**
 * @brief set the percentage for reporting #LOG_PROG messages.
 */
void Logging::report_percentage (const int percentage)
{
    if (percentage <= 0) {
        LOG_WARN << "Logging report percentage less than 1% not possible.";
        report_percentage_ = 1;
        return;
    }
    if (percentage > 100) {
        LOG_WARN << "Logging report percentage greater than 100% not meaningful.";
        report_percentage_ = 100;
        return;
    }
    report_percentage_ = percentage;
}

/**
 * @brief Return a string representation of a log level.
 */
std::string Logging::LevelToString(const LoggingLevel level)
{
    static const char* const buffer[] = {
        "NONE", "ERR ", "WARN", "INFO", "PROG", "DBG ", "DBG1", "DBG2", "DBG3", "DBG4"
    };
    return buffer[level];
}

/**
 * @brief Add stdout as output stream to which log messages are written.
 */
void Logging::LogToStdout ()
{
    // check whether stdout was already added.
    for (std::ostream* os : ostreams_) {
        if (os == &std::cout) {
            return;
        }
    }

    // if not, add it as output stream.
    ostreams_.push_back (&std::cout);
}

/**
 * @brief Add an output stream to which log messages are written.
 */
void Logging::LogToStream (std::ostream& os)
{
    ostreams_.push_back (&os);
}

/**
 * @brief Add an output file to which log messages are written.
 *
 * This creates a stream to the file.
 */
void Logging::LogToFile (const std::string& fn)
{
    // TODO the log file stream is never deleted. this is not a big leak,
    // as commonly only one file is used for logging, but still is a smell.
    std::ofstream* file = new std::ofstream();
    file->open (fn, std::ios::out | std::ios::app);
    if (file->is_open()) {
        ostreams_.push_back (file);
    } else {
        LOG_WARN << "Cannot open logging file " << fn;
    }
}

// =============================================================================
//     Destructor (does the actual work)
// =============================================================================

// TODO the output of the log is not thread safe
/**
 * @brief Destructor that is invoked at the end of each log line and does the actual
 * output.
 */
Logging::~Logging()
{
    // build the details for the log message into a buffer
    clock_t now_clock = clock();
    std::ostringstream det_buff;
    det_buff.str("");
    if (details_.count) {
        det_buff.fill('0');
        det_buff.width(4);
        det_buff << count_ << " ";
    }
    if (details_.date) {
        det_buff << CurrentDate() << " ";
    }
    if (details_.time) {
        det_buff << CurrentTime() << " ";
    }
    if (details_.runtime) {
        det_buff << std::fixed
                 << std::setprecision(6)
                 << double(now_clock) / CLOCKS_PER_SEC
                 << " ";
    }
    if (details_.rundiff) {
        double val = 0.0;
        if (last_clock_ > 0) {
            val = (double) (now_clock - last_clock_) / CLOCKS_PER_SEC;
        }
        det_buff << std::fixed
                 << std::setprecision(6)
                 << val
                 << " ";
        last_clock_ = now_clock;
    }
    if (details_.file) {
        det_buff << file_ << (details_.line ? "" : " ");
    }
    if (details_.line) {
        det_buff << ":" << line_ << " ";
    }
    if (details_.function) {
        det_buff << "(" << function_ << ") ";
    }
    if (details_.level) {
        det_buff << LevelToString(level_) << " ";
    }

    // add spaces for nested debug levels
    if (level_ > kDebug) {
        for (int i = 0; i < level_ - kDebug; i++) {
            det_buff << debug_indent;
        }
    }

    // make multi line log messages align to the length of the detail header,
    // and trim trailing whitespace, as we only want one newline at the end
    std::string msg = det_buff.str();
    if (msg.length() > 0) {
        msg += StringReplaceAll(
            buff_.str(), "\n", "\n" + std::string(msg.length(), ' ')
        );
    } else {
        msg += buff_.str();
    }
    msg = StringTrimRight(msg);

    // output the message to every stream, thread safe!
#   ifdef PTHREADS
    log_mutex.lock();
#   endif
    for (std::ostream* out : ostreams_) {
        (*out) << msg << std::endl << std::flush;
    }
#   ifdef PTHREADS
    log_mutex.unlock();
#   endif

    // inc log message counter
    count_++;
}

// =============================================================================
//     Singleton accessors
// =============================================================================

/**
 * @brief Getter for the singleton instance of log, is called by the standard macros.
 *
 * It returns the string stream buffer used to capture the log messages.
 */
std::ostringstream& Logging::Get(
    const std::string& file, const int line, const std::string& function,
    const LoggingLevel level
)
{
    return Get(file, line, function, level, details);
}

/**
 * @brief Getter for the singleton instance of log, is called by special macros
 * that change the details of the log message.
 *
 * It stores some relevant information and returns the string stream buffer
 * used to capture the log messages.
 */
std::ostringstream& Logging::Get(
    const std::string& file, const int line, const std::string& function,
    const LoggingLevel level, const LoggingDetails dets
)
{
    // save the information given when called from the macros
    file_     = file;
    line_     = line;
    function_ = function;
    level_    = level;
    details_  = dets;
    buff_.str("");
    return buff_;
}

} // namespace genesis
