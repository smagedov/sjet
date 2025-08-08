#ifndef SJET_CLUSTERINGSEQUENCE_HH_
#define SJET_CLUSTERINGSEQUENCE_HH_

#include <cmath>
#include <cassert>
#include <vector>
#include <utility>
#include <set>
#include <limits>

#include "sjet/Cluster.hh"

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

        DistCalc distCalc_;
        std::vector<cluster_type> clustHist_;
        DistSet distSet_;
        double maxClusDist_;
        int nParticles_;
        int nClusters_;
    };
}

#endif // SJET_CLUSTERINGSEQUENCE_HH_
