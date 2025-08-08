#ifndef SJET_CLUSTERINGSEQUENCE_HH_
#define SJET_CLUSTERINGSEQUENCE_HH_

#include <cmath>
#include <cassert>
#include <vector>
#include <utility>
#include <set>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "sjet/Cluster.hh"
#include "sjet/AbsNodeVisitor.hh"

namespace sjet {
    template <class DistCalc, class Particle>
    class ClusteringSequence {
    public:
        typedef DistCalc calc_type;
        typedef Particle particle_type;
        typedef Cluster<Particle> cluster_type;

        ClusteringSequence(const DistCalc& calc)
            : distCalc_(calc) {reset();}

        inline int nParticles() const {return nParticles_;}
        inline int nClusters() const {return nClusters_;}
        inline const std::vector<cluster_type>& clustHist() const
            {return clustHist_;}
        inline double maxCusteredDistance() const {return maxClusDist_;}

        void init(const std::vector<Particle>& initialParts) {
            reset();
            nParticles_ = initialParts.size();
            nClusters_ = nParticles_;
            if (nParticles_ > 0) {
                clustHist_.reserve(2*nParticles_-1);
                for (int i = 0; i < nParticles_; ++i) {
                    clustHist_.emplace_back(initialParts[i]);
                }

                for (int i = 1; i < nParticles_; ++i) {
                    for (int j = 0; j < i; ++j) {
                        const double dist = distCalc_(initialParts[j], initialParts[i]);
                        assert(dist >= 0.0);
                        const bool status = distSet_.insert(DistElem(dist, IndexPair(j, i))).second;
                        assert(status);
                    }
                }
            }
        }

        bool recomb() {
            if (distSet_.empty()) {
                if (nParticles_) {
                    assert(1 == nClusters_);
                }
                return false;
            }
            const DistSet::const_iterator it = distSet_.begin();
            const int newind = clustHist_.size();
            const int p1 = it->second.first;
            const int p2 = it->second.second;
            if (it->first > maxClusDist_)
                maxClusDist_ = it->first;
            clustHist_.emplace_back(clustHist_[p1], p1, clustHist_[p2],
                                    p2, it->first, maxClusDist_);
            clustHist_[p1].setDaughter(newind);
            clustHist_[p2].setDaughter(newind);
            updateDistTable(newind);
            tableClean();
            --nClusters_;
            return true;
        }

        inline double nextDistance() const {
            if (distSet_.empty())  {
                return -1.0;
            } else {
                return distSet_.begin()->first;
            }
        }

        inline double lastDistance() const {
            if (clustHist_.empty()) {
                return -1.0;
            } else {
                return clustHist_.back().dist();
            }
        }

        void run(const double maxDistance = std::numeric_limits<double>::max()) {
            while (nClusters_ > 1) {
                if (nextDistance() <= maxDistance) {
                    recomb();
                } else {
                    break;
                }
            }
        }

        // Return the vector of recombination distances
        std::vector<double> recombinationDistances() const
        {
            std::vector<double> distances;
            const int sz = clustHist_.size();
            if (sz > nParticles_)
            {
                distances.reserve(sz - nParticles_);
                for (int i=nParticles_; i<sz; ++i)
                    distances.push_back(clustHist_[i].dist());
            }
            return distances;
        }

        // Indices of clusters without daughters (final clusters)
        // in the complete clustering history. These are meaningful
        // if the clustering was run to some max distance and not
        // until only one cluster remains.
        std::vector<unsigned> clusterIndices() const
        {
            std::vector<unsigned> result;
            result.reserve(nClusters_);
            const unsigned lenHist = clustHist_.size();
            for (unsigned i=0; i<lenHist; ++i)
                if (!clustHist_[i].daughter())
                    result.push_back(i);
            return result;
        }

        // Clusters that we get at a particular distance (even if
        // the clustering was run until only one cluster remains)
        std::vector<unsigned> clusterIndices(const double dist) const
        {
            assert(dist > 0.0);

            std::vector<unsigned> result;
            result.reserve(nClusters_);
            const unsigned lenHist = clustHist_.size();
            if (lenHist)
            {
                if (dist < maxClusDist_)
                {
                    const unsigned maxClus = firstBigCluster(dist);
                    for (unsigned i=0; i<maxClus; ++i)
                    {
                        const unsigned dau = clustHist_[i].daughter();
                        if (!dau || dau >= maxClus)
                            result.push_back(i);
                    }
                }
                else
                {
                    // Did we run the clustering sequence at least
                    // until "dist"?
                    if (1U == nClusters_)
                    {
                        // We can assume that the clustering sequence
                        // was run to the end
                        result.push_back(lenHist - 1U);
                    }
                    else if (dist == maxClusDist_ || dist < nextDistance())
                    {
                        // The easy case: just collect clusters without
                        // daughters
                        for (unsigned i=0; i<lenHist; ++i)
                            if (!clustHist_[i].daughter())
                                result.push_back(i);
                    }
                    else
                    {
                        // No, we did not run clustering until "dist".
                        // In such a situation we can't determine the
                        // cluster set at this distance.
                        throw std::runtime_error(run_msg("clusterIndices", dist));
                    }
                }
            }
            return result;
        }

        // Cluster assignments of initial particles in the
        // complete clustering history
        std::vector<unsigned> clusterAssignments() const
        {
            std::vector<unsigned> result(nParticles_);
            const int lenHist = clustHist_.size();
            for (int i=0; i<nParticles_; ++i)
                result[i] = topDau(i, lenHist);
            return result;
        }

        // Cluster assignments of initial particles
        // that we get at a particular distance (even if
        // the clustering was run until only one cluster remains)
        std::vector<unsigned> clusterAssignments(const double dist) const
        {
            std::vector<unsigned> result(nParticles_);
            const int lenHist = clustHist_.size();
            if (lenHist)
            {
                if (dist < maxClusDist_)
                {
                    const int maxClus = firstBigCluster(dist);
                    for (int i=0; i<nParticles_; ++i)
                        result[i] = topDau(i, maxClus);
                }
                else
                {
                    // Did we run the clustering sequence at least
                    // until "dist"?
                    if (1U == nClusters_)
                    {
                        // We can assume that the clustering sequence
                        // was run to the end
                        std::fill(result.begin(), result.end(), lenHist - 1U);
                    }
                    else if (dist == maxClusDist_ || dist < nextDistance())
                    {
                        for (int i=0; i<nParticles_; ++i)
                            result[i] = topDau(i, lenHist);
                    }
                    else
                    {
                        // No, we did not run clustering until "dist".
                        // In such a situation we can't determine the
                        // clustering assignments at this distance.
                        throw std::runtime_error(run_msg("clusterAssignments", dist));
                    }
                }
            }
            return result;
        }

        // Visit all particles (skipping intermediate nodes)
        // that are recombined into the current node
        void visitParentParticles(AbsNodeVisitor<Particle>& visitor,
                                  const int node) const
        {
            const cluster_type& clus = clustHist_.at(node);
            if (clus.hasParents())
            {
                visitParentParticles(visitor, clus.parent1());
                visitParentParticles(visitor, clus.parent2());
            }
            else
                visitor.visit(static_cast<unsigned>(node), clus.p());
        }

        // Visit the given node and all of its daughters,
        // up to some recombination distance cutoff
        void visitNodeAndDaus(AbsNodeVisitor<Particle>& visitor,
                              const int node, const double upToDist) const
        {
            const cluster_type& clus = clustHist_.at(node);
            if (upToDist < clus.distance())
                return;
            const int maxClus = upToDist < maxClusDist_ ?
                                firstBigCluster(upToDist) :
                                static_cast<int>(clustHist_.size());
            visitNodeAndDaus2(visitor, node, maxClus);
        }

        // Visit the ancestors of this cluster, choosing each time
        // the parent which carries the "larger" particle.
        // Instance of class Comp should work like operator "less"
        // for particles.
        template<class Comp>
        void visitLargerAncestors(AbsNodeVisitor<Particle>& visitor,
                                  const int node,
                                  const Comp& parentComparator) const
        {
            const cluster_type& clus = clustHist_.at(node);
            if (clus.hasParents())
            {
                const int par1 = clus.parent1();
                const int par2 = clus.parent2();
                const Particle& p1 = clustHist_[par1].p();
                const Particle& p2 = clustHist_[par2].p();
                if (parentComparator(p1, p2)) // p1 compares less than p2
                {
                    visitor.visit(static_cast<unsigned>(par2), p2);
                    visitLargerAncestors(visitor, par2, parentComparator);
                }
                else
                {
                    visitor.visit(static_cast<unsigned>(par1), p1);
                    visitLargerAncestors(visitor, par1, parentComparator);
                }
            }
        }

    private:
        typedef std::pair<int, int> IndexPair;
        typedef std::pair<double, IndexPair> DistElem;
        typedef std::set<DistElem> DistSet;

        void reset() {
            clustHist_.clear();
            distSet_.clear();
            maxClusDist_ = -1.0;
            nParticles_ = 0;
            nClusters_ = 0;
        }

        void updateDistTable(const int newind) {
            const Particle& newpart = clustHist_[newind].p();
            for (int i = 0; i < newind; ++i) {
                if (!clustHist_[i].daughter()) {
                    const double dist = distCalc_(clustHist_[i].p(), newpart);
                    assert(dist >= 0.0);
                    const bool status = distSet_.insert(DistElem(dist, IndexPair(i, newind))).second;
                    assert(status);
                }
            }
        }

        void tableClean() {
            while (!distSet_.empty()) {
                DistSet::iterator it = distSet_.begin();
                if (clustHist_[it->second.first].daughter() ||
                    clustHist_[it->second.second].daughter()) {
                    distSet_.erase(it);
                } else {
                    break;
                }
            }
        }

        int topDau(const int particleInd, const int maxClus) const
        {
            assert(particleInd < maxClus);
            const int dau = clustHist_[particleInd].daughter();
            if (!dau || dau >= maxClus)
                return particleInd;
            else
                return topDau(dau, maxClus);
        }

        std::string run_msg(const std::string& where, const double dist) const
        {
            std::ostringstream os;
            os.precision(16);
            os << "In sjet::ClusteringSequence::" << where
               << " : distance argument " << dist << " is too large."
               << " In order for this argument to work, run the"
               << " clustering sequence with increased distance cutoff.";
            return os.str();
        }

        // Find the first cluster whose recombination
        // distance is larger than the argument
        unsigned firstBigCluster(const double dist) const
        {
            assert(dist > 0.0);
            const unsigned lenHist = clustHist_.size();
            if (dist < maxClusDist_)
                for (unsigned i=nParticles_; i<lenHist; ++i)
                    if (clustHist_[i].dist() > dist)
                        return i;
            return lenHist;
        }

        void visitNodeAndDaus2(AbsNodeVisitor<Particle>& visitor,
                               const int node, const int maxClus) const
        {
            const cluster_type& clus = clustHist_[node];
            visitor.visit(static_cast<unsigned>(node), clus.p());
            const int dau = clus.daughter();
            if (dau && dau < maxClus)
                visitNodeAndDaus2(visitor, dau, maxClus);
        }

        DistCalc distCalc_;
        std::vector<cluster_type> clustHist_;
        DistSet distSet_;
        double maxClusDist_;
        int nParticles_;
        int nClusters_;
    };
}

#endif // SJET_CLUSTERINGSEQUENCE_HH_
