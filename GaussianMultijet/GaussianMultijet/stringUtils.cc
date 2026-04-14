#include <sstream>
#include <ctime>
#include <cstdio>

#include "stringUtils.hh"

std::string appendAnInteger(const std::string& prefix, const long l,
                            const char sep)
{
    std::ostringstream os;
    os << prefix;
    if (sep != '\0')
        os << sep;
    os << l;
    return os.str();
}

std::string concatStringsAndInts(const std::string& prefix,
                                 const long i1, const long i2,
                                 const std::string& postfix,
                                 const char sep)
{
    std::ostringstream os;
    os << prefix;
    if (sep != '\0')
        os << sep;
    os << i1;
    if (sep != '\0')
        os << sep;
    os << i2;
    os << postfix;
    return os.str();
}

std::string timestamp()
{
    struct tm *current;
    time_t now;
    
    time(&now);
    current = localtime(&now);

    char buf[16];
    sprintf(buf, "%02i:%02i:%02i", current->tm_hour,
            current->tm_min, current->tm_sec);
    return std::string(buf);
}
