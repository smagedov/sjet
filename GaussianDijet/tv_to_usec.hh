#ifndef TV_TO_USEC_HH_
#define TV_TO_USEC_HH_

#include <sys/time.h>

inline unsigned long long tv_to_usec(const struct timeval* tv)
{
    unsigned long long usec = tv->tv_usec;
    usec += tv->tv_sec*1000000ULL;
    return usec;
}

#endif // TV_TO_USEC_HH_
