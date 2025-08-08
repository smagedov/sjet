#ifndef TRANSVERSEMMD_HH_
#define TRANSVERSEMMD_HH_

#include <set>
#include <map>
#include <vector>
#include <algorithm>

#include "MatrixMapper.hh"
#include "permutation.hh"

template <class Particle>
double transverseMMD(const std::vector<Particle>& particles,
                     const std::vector<unsigned>& clustAssignments1,
                     const std::vector<unsigned>& clustAssignments2)
{
    const unsigned nParticles = particles.size();
    assert(nParticles);
    assert(clustAssignments1.size() == nParticles);
    assert(clustAssignments2.size() == nParticles);

    const std::set<unsigned> idSet1(clustAssignments1.begin(), clustAssignments1.end());
    const unsigned nClus1 = idSet1.size();
    std::map<unsigned,unsigned> lookup1;
    {
        unsigned clusNum = 0;
        for (std::set<unsigned>::const_iterator it = idSet1.begin(); it != idSet1.end(); ++it)
            lookup1[*it] = clusNum++;
    }

    const std::set<unsigned> idSet2(clustAssignments2.begin(), clustAssignments2.end());
    const unsigned nClus2 = idSet2.size();
    std::map<unsigned,unsigned> lookup2;
    {
        unsigned clusNum = 0;
        for (std::set<unsigned>::const_iterator it = idSet2.begin(); it != idSet2.end(); ++it)
            lookup2[*it] = clusNum++;
    }

    const unsigned nClusMax = std::max(nClus1, nClus2);
    const MatrixMapper m(nClusMax, nClusMax);
    std::vector<long double> mat(m.size(), 0.0L);
    long double ptSum = 0.0L;
    for (unsigned i=0; i<nParticles; ++i)
    {
        const double pt = particles[i].pt();
        ptSum += pt;
        mat[m(lookup1[clustAssignments1[i]], lookup2[clustAssignments2[i]])] += pt;
    }

    // Brute force algorithm which cycles over all
    // possible permutations (there are nClusMax! of them).
    // Of course, this will be slow for large nClusMax.
    std::vector<unsigned> perm(nClusMax);
    long double maxPtSum = 0.0L;
    const unsigned long nPerms = factorial(nClusMax);
    for (unsigned long i=0; i<nPerms; ++i)
    {
        orderedPermutation(i, &perm[0], nClusMax);
        long double permSum = 0.0L;
        for (unsigned row=0; row<nClusMax; ++row)
            permSum += mat[m(row, perm[row])];
        if (permSum > maxPtSum)
            maxPtSum = permSum;
    }

    return 1.0L - maxPtSum/ptSum;
}

#endif // TRANSVERSEMMD_HH_
