#ifndef SM_CLUSTERINGSEQUENCE_HH_
#define SM_CLUSTERINGSEQUENCE_HH_

#include <cmath>
#include <cassert>
#include <vector>
#include <utility>
#include <set>
#include <limits>

namespace sjet {

	template <class Particle>
	class Cluster {
		public:
			Cluster(const Particle& originpart) 
				: p_(originpart), dist_(-1.0), parent1_(-1), parent2_(-1), daughter_(0) {}

			Cluster(const Particle& p1, const int index1, 
				const Particle& p2, const int index2, const double distance)
				: p_(p1 + p2), dist_(distance), parent1_(index1), 
				parent2_(index2), daughter_(0) {
				assert(parent1_ >= 0);
				assert(parent2_ >= 0);
				assert(dist_ >= 0.0);
				}
			
			inline double dist() const {return dist_;}
			inline const Particle& p() const {return p_;}
			inline int parent1() const {return parent1_;}
			inline int parent2() const {return parent2_;}
			inline int daughter() const {return daughter_;}

			inline void setDaughter(const int d) {
				assert(d > 0);
				daughter_ = d;
			};


		private:
			Particle p_;
			double dist_;
			int parent1_;
			int parent2_;
			int daughter_;

	};

	template <class DistCalc, class Particle>
	class ClusteringSequence {
		public:
			ClusteringSequence(const DistCalc& calc)
			       : distCalc_(calc) {
				       reset();
			       }

			inline int nClusters() const {return nClusters_;}
			inline std::vector<Cluster<Particle> > clustHist() const {return clustHist_;}

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
				clustHist_.emplace_back(clustHist_[p1].p(), p1, clustHist_[p2].p(), p2, it->first);
				clustHist_[p1].setDaughter(newind);
				clustHist_[p2].setDaughter(newind);
				updateDistTable(newind);

//				std::cout << nClusters_ << std::endl;
				
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
					if (nextDistance() < maxDistance) {
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
					const int p1 = it->second.first;
                                	const int p2 = it->second.second;
					if (clustHist_[p1].daughter() || clustHist_[p2].daughter()) {
						distSet_.erase(it);
					} else {
						break;
					}
				}
			}

			DistCalc distCalc_;
			std::vector<Cluster<Particle> > clustHist_;
			DistSet distSet_;
			int nParticles_;
			int nClusters_;
	};

}

#endif // SM_CLUSTERINGSEQUENCE_HH_
