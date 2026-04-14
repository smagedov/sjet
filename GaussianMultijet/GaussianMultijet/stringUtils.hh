#ifndef STRINGUTILS_HH_
#define STRINGUTILS_HH_

#include <string>

// A time stamping routine for simple printouts. Uses hh:mm:ss format.
std::string timestamp();

// The default corresponds to not having a separator character
// between the prefix and the appended integer
std::string appendAnInteger(const std::string& prefix, long l,
                            char separator = '\0');

// This will concatenate prefix, separator, i1, separator, i2, postfix
std::string concatStringsAndInts(const std::string& prefix,
                                 long i1, long i2,
                                 const std::string& postfix,
                                 char separator = '_');

#endif // STRING_UTILS_HH_
