#include <numeric>

#include "genClusterAssignments.hh"

std::vector<unsigned> genClusterAssignments(
    const std::vector<unsigned>& particlesPerCluster)
{
    const unsigned nParticles = std::accumulate(
        particlesPerCluster.begin(), particlesPerCluster.end(), 0U);
    std::vector<unsigned> result(nParticles);
    const unsigned nClus = particlesPerCluster.size();
    unsigned imin = 0;
    for (unsigned iclus = 0; iclus < nClus; ++iclus)
    {
        const unsigned imax = imin + particlesPerCluster[iclus];
        for (unsigned i=imin; i<imax; ++i)
            result[i] = iclus;
        imin = imax;
    }
    return result;
}
