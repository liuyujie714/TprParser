#ifndef TPR_EXCEPTION_H
#define TPR_EXCEPTION_H

#include <stdexcept>
#include <string>

class _Error : public std::runtime_error
{
public:
    _Error(const std::string& msg,
        const char* fpath,
        int line,
        const char* funcname)
        : runtime_error(msg)
    {
        msg_ = "Error: " + msg + "\n\nSource File:\t" + fpath + ", Line: " + std::to_string(line) +
            "\nFunction:\t" + funcname;
    }

    virtual const char* what() const noexcept override
    {
        return msg_.c_str();
    }

private:
    std::string			msg_;
};

#ifdef _MSC_VER
//! throw error
#define THROW_TPR_EXCEPTION(msg) throw _Error(msg, __FILE__, __LINE__, __FUNCSIG__)
#elif defined(__PRETTY_FUNCTION__)
//! throw error
#define THROW_TPR_EXCEPTION(msg) _Error(msg, __FILE__, __LINE__, __PRETTY_FUNCTION__)
#else
#define THROW_TPR_EXCEPTION(msg) throw _Error(msg, __FILE__, __LINE__, "Unknown")
#endif

#endif // !TPR_EXCEPTION_H
