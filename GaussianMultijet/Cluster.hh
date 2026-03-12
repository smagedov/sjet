#ifndef SC_CLUSTER_HH_
#define SC_CLUSTER_HH_

#include <cassert>

namespace sc {
    template<class Particle>
    class Cluster
    {
    public:
        typedef Particle particle_class;

        // Construct a cluster consisting of one particle
        inline explicit Cluster(const Particle& p)
            : particle_(p), d_(-1.0), maxd_(-1.0),
              scalarPtSum_(p.pt()),
              p1ind_(-1), p2ind_(-1), dau_(0) {}

        // Construct a cluster consisting of two particles
        inline Cluster(const Cluster& c1, const int ind1,
                       const Cluster& c2, const int ind2,
                       const double dist, const double maxdist)
            : particle_(c1.particle() + c2.particle()),
              d_(dist), maxd_(maxdist),
              scalarPtSum_(c1.scalarPtSum() + c2.scalarPtSum()),
              p1ind_(ind1), p2ind_(ind2), dau_(0)
        {
            assert(p1ind_ >= 0);
            assert(p2ind_ >= 0);
            assert(d_ >= 0.0);
            assert(maxd_ >= 0.0);
        }

        // Modifiers
        inline void setDau(const int dau)
        {
            assert(dau > 0);
            dau_ = dau;
        }

        // Inspectors
        inline const Particle& particle() const {return particle_;}
        inline double distance() const {return d_;}
        inline double maxDistanceSoFar() const {return maxd_;}
        inline double scalarPtSum() const {return scalarPtSum_;}
        inline int parent1() const {return p1ind_;}
        inline int parent2() const {return p2ind_;}
        inline int dau() const {return dau_;}
        inline bool hasParents() const
        {
            if (p1ind_ >= 0)
            {
                assert(p2ind_ >= 0);
                return true;
            }
            else
            {
                assert(p2ind_ < 0);
                return false;
            }
        }
        inline bool hasDau() const {return dau_ > 0;}

        inline bool operator==(const Cluster& r) const
        {
            return particle_ == r.particle_ &&
                   d_ == r.d_ &&
                   maxd_ == r.maxd_ &&
                   scalarPtSum_ == r.scalarPtSum_ &&
                   p1ind_ == r.p1ind_ &&
                   p2ind_ == r.p2ind_ &&
                   dau_ == r.dau_;
        }

    private:
        Particle particle_;
        double d_;    // Distance at recombination, -1 if this is an original
                      // particle not created in the recombination process
        double maxd_; // Maximum distance observed so far in the recombination
                      // sequence, -1 if this is an original particle
        double scalarPtSum_;
        int p1ind_;   // Parent 1 index, -1 if none
        int p2ind_;   // Parent 2 index, -1 if none
        int dau_;     // Daughter index,  0 if none
    };
}

#endif // SC_CLUSTER_HH_
