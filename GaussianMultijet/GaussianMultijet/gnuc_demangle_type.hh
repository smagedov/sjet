#ifndef GNUC_DEMANGLE_TYPE_HH_
#define GNUC_DEMANGLE_TYPE_HH_

#include <string>

#ifdef __GNUC__

#include <cstdlib>
#include <typeinfo>
#include <cxxabi.h>

template <typename T>
inline std::string gnuc_demangle_type(const T& obj)
{
    int status;
    const std::type_info &ti = typeid(obj);
    char *realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
    std::string s(realname);
    free(realname);
    return s;
}

#else

template <typename T>
inline std::string gnuc_demangle_type(const T&)
{
    return "Error: unable to demangle type name for this compiler";
}

#endif // __GNUC__

#endif // GNUC_DEMANGLE_TYPE_HH_
