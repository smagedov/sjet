#ifndef DIJETEVENT_HH_
#define DIJETEVENT_HH_

#include <limits>
#include <stdexcept>

#include "sjet/ClusteringSequence.hh"
#include "sjet/DistanceCalculator.hh"

#include "DijetComponents.hh"

// When you add new items to the event, don't forget to add the
// corresponding "Ready" flag and set it to "false" inside the
// "clear" function.
template <class Particle>
struct DijetEvent
{
    typedef Particle particle_type;
    typedef sjet::ClusteringSequence<sjet::DistanceCalculator, Particle> clust_seq_type;

    inline DijetEvent()
        : diffusionSequence((sjet::DistanceCalculator())) {clear();}

    // Expected direction of the probe jet
    double expectedProbeEta;
    double expectedProbePhi;

    // Generated event contents
    std::pair<std::vector<Particle>, std::vector<unsigned> > genEvent;
    bool genEventReady;

    // 4-momenta of the generator jets
    std::vector<Particle> genJets;
    bool genJetsReady;

    // Clustering sequence for the diffusion distance
    clust_seq_type diffusionSequence;
    bool diffusionSequenceReady;

    // 4-momenta of the jets obtained with the diffusion distance
    // clustering and defined by a simple constant distance cutoff.
    // These 4-momenta will be ordered by pt in the decreasing order.
    std::vector<Particle> simpleDiffusionJets;

    // Corresponding nodes in the clustering history
    std::vector<unsigned> simpleDiffusionNodes;
    double simpleDiffusionJetsDistCutoff;
    bool simpleDiffusionJetsReady;

    // pt-weighted Rand and Jaccard clustering distances
    // as well as MMD (see the paper by Cowden and Volobouev)
    double simpleDiffusionRandDist;
    double simpleDiffusionJaccardDist;
    double simpleDiffusionMMDDist;
    bool simpleDiffusionClusDistReady;

    // This function should be called every time
    // before generating and processing a new event
    inline void clear()
    {
        number_ = invalidEventNumber_;
        genEventReady = false;
        genJetsReady = false;
        diffusionSequenceReady = false;
        simpleDiffusionJetsReady = false;
        simpleDiffusionClusDistReady = false;
    }

    // Helper functions calculating some dependent quantities
    // that may later be of interest. For simplicity of saving
    // the results, these functions should take no arguments.

    // Total multiplicity of the generated particles
    // (after the Pt cut, if any)
    inline unsigned genTotalMult() const
    {
        assert(genEventReady);
        return genEvent.first.size();
    }

    // Multiplicity of the generated pileup (if any)
    inline unsigned genPileupMult() const
    {
        assert(genEventReady);
        return genEvent.second.at(DC_PILEUP);
    }

    // Scalar Pt sum of the generated pileup particles
    inline double genPileupScalarPtSum() const
    {
        assert(genEventReady);
        const std::pair<unsigned, unsigned> minmax =
            dijetComponentRange(genEvent.second, DC_PILEUP);
        long double ptSum = 0.0L;
        for (unsigned i=minmax.first; i<minmax.second; ++i)
            ptSum += genEvent.first[i].pt();
        return ptSum;
    }

    // Standard functions -- see MinimalEvent.hh
    // Handle event numbers
    inline void setNumber(const unsigned long n)
    {
        if (n == invalidEventNumber_)
            throw std::invalid_argument("In DijetEvent::setNumber : "
                                        "invalid event number");
        number_ = n;
    }
    inline unsigned long number() const {return number_;}

    // Check if the event number has been set
    inline bool isNumberValid() const
        {return number_ != invalidEventNumber_;}

    // Bump up the version number if you change
    // event contents
    static inline unsigned version() {return 1;}

private:
    static constexpr unsigned long invalidEventNumber_ =
        std::numeric_limits<unsigned long>::max();
    unsigned long number_;
};

#endif // DIJETEVENT_HH_
