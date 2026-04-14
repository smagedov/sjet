#ifndef MATCHPAIRTOPAIR_HH_
#define MATCHPAIRTOPAIR_HH_

#include <utility>
#include <vector>
#include <cassert>
#include <limits>

//
// Find the best two elements among "candidates" that have the
// smallest summary distance to obj0 and obj1. The vector of
// candidates must have at least two elements. In the returned pair,
// the first element is the index in the vector of candidates
// for obj0 match, and the second is the index for obj1 match.
//
// The distance functor should have a method which looks like
// double operator()(const T& left, const T& right) const;
//
template<typename T, class DistanceFunctor>
std::pair<unsigned,unsigned> matchPairToPair(
    const T& obj0, const T& obj1, const std::vector<T>& candidates,
    const DistanceFunctor& distCalc)
{
    const unsigned nCand = candidates.size();
    assert(nCand > 1U);
    std::vector<double> buf(nCand*2U);
    double* dist0 = &buf[0];
    double* dist1 = dist0 + nCand;

    for (unsigned i=0; i<nCand; ++i)
    {
        dist0[i] = distCalc(obj0, candidates[i]);
        dist1[i] = distCalc(obj1, candidates[i]);
    }

    // Go over all possible pairs in "candidates" and
    // choose the pair with the smallest total distance
    std::pair<unsigned,unsigned> match(nCand, nCand);
    double bestDistance = std::numeric_limits<double>::max();

    for (unsigned cand0=0; cand0<nCand; ++cand0)
        for (unsigned cand1=0; cand1<nCand; ++cand1)
            if (cand1 != cand0)
            {
                const double dist = dist0[cand0] + dist1[cand1];
                if (dist < bestDistance)
                {
                    bestDistance = dist;
                    match.first = cand0;
                    match.second = cand1;
                }
            }

    return match;
}

#endif // MATCHPAIRTOPAIR_HH_
