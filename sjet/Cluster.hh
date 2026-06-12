#ifndef SJET_CLUSTER_HH_
#define SJET_CLUSTER_HH_

#include <cassert>
#include "objectPt.hh"

namespace sjet {
    template <class DistCalc, class SeqParticle>
    class ClusteringSequence;

    template <class Particle>
    class Cluster {
    public:
        typedef Particle particle_class;

        inline explicit Cluster(const Particle& originpart)
            : p_(originpart), dist_(-1.0), maxd_(-1.0),
              scalarPtSum_(objectPt(originpart)),
              parent1_(-1), parent2_(-1), daughter_(0), prob_(1), logprob_(0) {}

        // Copy-like constructor from another cluster type
        template<class OldParticle>
        inline explicit Cluster(const Cluster<OldParticle>& c, const Particle& newpart)
            : p_(newpart), dist_(c.dist_), maxd_(c.maxd_),
              scalarPtSum_(c.scalarPtSum_),
              parent1_(c.parent1_), parent2_(c.parent2_), daughter_(c.daughter_), 
	      prob_(c.prob_), logprob_(c.prob_) {}

        inline Cluster(const Cluster& c1, const int index1, 
                       const Cluster& c2, const int index2,
                       const double distance, const double maxdist)
            : p_(c1.p() + c2.p()), dist_(distance), maxd_(maxdist),
              scalarPtSum_(c1.scalarPtSum() + c2.scalarPtSum()),
              parent1_(index1), parent2_(index2), daughter_(0), 
	      prob_(), logprob_(c1.logprob() + c2.logprob()) {
            assert(parent1_ >= 0);
            assert(parent2_ >= 0);
            assert(dist_ >= 0.0);
            assert(maxd_ >= 0.0);
        }

        // Modifiers
        inline void setDaughter(const int d) {
            assert(d > 0);
            daughter_ = d;
        };

        // Inspectors
        inline const Particle& p() const {return p_;}
        inline double dist() const {return dist_;}
        inline double maxDistanceSoFar() const {return maxd_;}
        inline double scalarPtSum() const {return scalarPtSum_;}
        inline int parent1() const {return parent1_;}
        inline int parent2() const {return parent2_;}
        inline int daughter() const {return daughter_;}
        inline double prob() const {return prob_;}
        inline double logprob() const {return logprob_;}
        inline bool hasParents() const
        {
            if (parent1_ >= 0)
            {
                assert(parent2_ >= 0);
                return true;
            }
            else
            {
                assert(parent2_ < 0);
                return false;
            }
        }
        inline bool hasDau() const {return daughter_ > 0;}

        inline bool operator==(const Cluster& r) const
        {
            return p_ == r.p_ &&
                   dist_ == r.dist_ &&
                   maxd_ == r.maxd_ &&
                   scalarPtSum_ == r.scalarPtSum_ &&
                   parent1_ == r.parent1_ &&
                   parent2_ == r.parent2_ &&
                   daughter_ == r.daughter_ &&
		   prob_ == r.prob_ &&
		   logprob_ == r.logprob_;
        }

    private:
        Particle p_;
        double dist_;   // Distance at recombination, -1 if this is an original
                        // particle not created in the recombination process
        double maxd_;   // Maximum distance observed so far in the recombination
                        // sequence, -1 if this is an original particle
        double scalarPtSum_; // Sum of scalar pt values of all particles
                             // in the cluster
        double prob_;   // Probability
        double logprob_;// Log of the probability
        int parent1_;   // Parent 1 index, -1 if none
        int parent2_;   // Parent 2 index, -1 if none
        int daughter_;  // Daughter index,  0 if none
	
	template <class DistCalc, class SeqParticle>
	friend class ClusteringSequence;
    };
}

#endif // SJET_CLUSTER_HH_
