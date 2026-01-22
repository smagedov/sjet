#ifndef FILTERBYPT_HH_
#define FILTERBYPT_HH_

#include <utility>
#include <vector>
#include <cassert>

template <class Particle>
void filterByPt(
    const std::pair<std::vector<Particle>, std::vector<unsigned> >& from,
    std::pair<std::vector<Particle>, std::vector<unsigned> >* to,
    const double ptMin)
{
    assert(to);

    if (ptMin <= 0.0)
        *to = from;
    else
    {
        to->first.clear();
        to->second.clear();
        to->second.reserve(from.second.size());

        unsigned imin = 0;
        for (unsigned blockLength : from.second)
        {
            unsigned count = 0;
            const unsigned imax = imin + blockLength;
            for (unsigned i=imin; i<imax; ++i)
                if (from.first.at(i).pt() >= ptMin)
                {
                    to->first.push_back(from.first[i]);
                    ++count;
                }
            imin = imax;
            to->second.push_back(count);
        }
    }
}

#endif // FILTERBYPT_HH_
