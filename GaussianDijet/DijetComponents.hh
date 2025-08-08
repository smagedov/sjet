#ifndef DIJETCOMPONENTS_HH_
#define DIJETCOMPONENTS_HH_

#include <cassert>
#include <vector>
#include <utility>

// Enum for the three components of the pseudo dijet event
enum DijetComponents {
    DC_REFERENCE_JET = 0,
    DC_PROBE_JET,
    DC_PILEUP,
    DC_N_EVENT_COMPONENTS
};

inline std::pair<unsigned, unsigned> dijetComponentRange(
    const std::vector<unsigned>& componentCounts, const unsigned which)
{
    assert(componentCounts.size() == DC_N_EVENT_COMPONENTS);
    assert(which < DC_N_EVENT_COMPONENTS);

    unsigned imin = 0, imax = 0;
    for (unsigned i=0; i<DC_N_EVENT_COMPONENTS; ++i)
    {
        imax = imin + componentCounts[i];
        if (which == i)
            break;
        imin = imax;
    }
    return std::pair<unsigned, unsigned>(imin, imax);
}

#endif // DIJETCOMPONENTS_HH_
