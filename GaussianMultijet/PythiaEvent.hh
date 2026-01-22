#ifndef MINIMALEVENT_HH_
#define MINIMALEVENT_HH_

#include <limits>
#include <stdexcept>
#include "Pythia8/Pythia.h"

#include "sjet/ClusteringSequence.hh"
#include "sjet/DistanceCalculator.hh"

#include "MultijetComponents.hh"

// Event structures playing well with the analysis framework
// should have the methods described below:
//
// Default constructor.
//
// Static method "version".
//
// Methods for dealing with event numbers: "setNumber", "number",
// "isNumberValid".
//
// Method "clear" for invalidating the event contents. It is expected
// that the event structure will be created only once before entering
// the event loop and then the "clear" and "setNumber" methods will
// be called inside the loop before everything else happens.
//
template <class Particle>
struct PythiaEvent
{
    typedef Particle particle_type;
    typedef sjet::ClusteringSequence<sjet::DistanceCalculator, Particle> clust_seq_type;
    
    inline PythiaEvent()
        : diffusionSequence((sjet::DistanceCalculator())) {clear();}


    // Data contained in this event (just an example)
    //double randomDatum;
    //bool randomDatumReady;
    Pythia8::Event * pythiaEvent;
    std::shared_ptr<Pythia8::Pythia> pythia;
    bool pythiaEventReady;
    std::vector<int> promptPartons;
    std::vector<std::vector<int>> genClusters;

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

    // Handle event numbers
    inline void setNumber(const unsigned long n)
    {
        if (n == invalidEventNumber_)
            throw std::invalid_argument("In MinimalEvent::setNumber : "
                                        "invalid event number");
        number_ = n;
    }
    inline unsigned long number() const {return number_;}

    // Check if the event number has been set
    inline bool isNumberValid() const
        {return number_ != invalidEventNumber_;}

    // Invalidate event number and various data sections
    void clear()
    {
	pythiaEvent = nullptr;
        number_ = invalidEventNumber_;
        pythiaEventReady = false;
        genJetsReady = false;
        genEventReady = false;
        diffusionSequenceReady = false;
        simpleDiffusionJetsReady = false;
        simpleDiffusionClusDistReady = false;


    }

    // Bump up the version number if you change
    // event contents
    static inline unsigned version() {return 1;}

private:
    static constexpr unsigned long invalidEventNumber_ =
        std::numeric_limits<unsigned long>::max();
    unsigned long number_;
};

#endif // MINIMALEVENT_HH_
