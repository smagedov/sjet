#ifndef TRANSVERSERANDJACCARD_HH_
#define TRANSVERSERANDJACCARD_HH_

#include <vector>
#include <cassert>
#include <utility>

// The first element of the returned pair is the pt-weighted
// Rand distance (that is, 1 - Rand index) and the second is
// the pt-weighted Jaccard distance (1 - Jaccard index).
// pt weighting makes these quantities infrared safe.
template <class Particle>
std::pair<double, double> transverseRandJaccard(
    const std::vector<Particle>& particles,
    const std::vector<unsigned>& clustAssignments1,
    const std::vector<unsigned>& clustAssignments2)
{
    long double sums[4] = {0.0L};

    const unsigned nParticles = particles.size();
    assert(nParticles > 1U);
    assert(clustAssignments1.size() == nParticles);
    assert(clustAssignments2.size() == nParticles);

    // Cycle over all particle pairs
    for (unsigned i=1; i<nParticles; ++i)
    {
        const double ptI = particles[i].pt();
        const unsigned cI1 = clustAssignments1[i];
        const unsigned cI2 = clustAssignments2[i];
            
        for (unsigned j=0; j<i; ++j)
        {
            const double weight = ptI*particles[j].pt();
            const bool same_in_1 = cI1 == clustAssignments1[j];
            const bool same_in_2 = cI2 == clustAssignments2[j];

            if (same_in_1 && same_in_2)
                sums[0] += weight;
            else if (same_in_1 && !same_in_2)
                sums[2] += weight;
            else if (!same_in_1 && same_in_2)
                sums[3] += weight;
            else
                sums[1] += weight;
        }
    }

    const long double jacDenom = sums[0] + sums[2] + sums[3];
    const long double total = jacDenom + sums[1];
    assert(total > 0.0L);

    const long double randIndex = (sums[0] + sums[1])/total;
    const long double jaccardIndex = jacDenom > 0.0L ? sums[0]/jacDenom : 0.0L;

    return std::pair<double, double>(1.0L - randIndex, 1.0L - jaccardIndex);
}

#endif // TRANSVERSERANDJACCARD_HH_
