#ifndef PARTICLESPECIALIZATIONS_HH_
#define PARTICLESPECIALIZATIONS_HH_

#include <stdexcept>
#include <vector>

#include "sjet/objectPt.hh"

template<class Particle>
std::vector<Particle> operator+(const std::vector<Particle>& l,
                                const std::vector<Particle>& r)
{
    const unsigned long sz = l.size();
    if (sz != r.size())
        throw std::invalid_argument("Incompatible vector sizes");
    std::vector<Particle> result;
    result.reserve(sz);
    for (unsigned long i=0; i<sz; ++i)
        result.push_back(l[i] + r[i]);
    return result;
}

namespace sjet {
    template<class Particle>
    double objectPt(const std::vector<Particle>& particles)
    {
        long double ptSum = 0.0;
        for (const auto& p : particles)
            ptSum += p.pt();
        return ptSum;
    }
}

#endif // PARTICLESPECIALIZATIONS_HH_
