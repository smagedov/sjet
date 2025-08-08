#ifndef RANDJACCARDMMDCALCULATOR_HH_
#define RANDJACCARDMMDCALCULATOR_HH_

#include "AbsFrameworkModule.hh"
#include "genClusterAssignments.hh"
#include "transverseRandJaccard.hh"
#include "transverseMMD.hh"
#include "DijetComponents.hh"

template <class Event>
struct RandJaccardMMDCalculator : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline RandJaccardMMDCalculator(const std::string& i_label,
                                    const double distanceCutoff)
        : Base(i_label), distanceCutoff_(distanceCutoff)
    {
        assert(distanceCutoff_ > 0.0);
    }

    inline virtual ~RandJaccardMMDCalculator() override {}

    inline virtual RandJaccardMMDCalculator* clone() const override
        {return new RandJaccardMMDCalculator(*this);}

    inline virtual bool process(Event& evt) override
    {
        // Make sure that we have the input data structures ready
        assert(evt.genEventReady);
        assert(evt.diffusionSequenceReady);
        assert(evt.simpleDiffusionJetsReady);

        // Make sure that we are not overwriting a part of the
        // event already made by some other module
        assert(!evt.simpleDiffusionClusDistReady);

        // Clustering assignments in the generated event
        const std::vector<unsigned>& genClusAssign =
            genClusterAssignments(evt.genEvent.second);

        // Clustering assignments in the diffusion clustering
        std::vector<unsigned> jetClusAssign =
            evt.diffusionSequence.clusterAssignments(distanceCutoff_);

        // Only use the diffusion clustering assignments to the
        // two leading pt jets. Assume everything else is pileup
        // or unclustered.
        const unsigned nDiffusionJets = evt.simpleDiffusionJets.size();
        if (nDiffusionJets > DC_N_EVENT_COMPONENTS)
        {
            const unsigned node0 = evt.simpleDiffusionNodes[0];
            const unsigned node1 = evt.simpleDiffusionNodes[1];
            const unsigned unusedNode = evt.diffusionSequence.clustHist().size();

            const unsigned sz = jetClusAssign.size();
            for (unsigned i=0; i<sz; ++i)
            {
                const unsigned a = jetClusAssign[i];
                if (a != node0 && a != node1)
                    jetClusAssign[i] = unusedNode;
            }
        }

        // Calculate the distances
        const std::pair<double, double>& distances =
            transverseRandJaccard(evt.genEvent.first, genClusAssign, jetClusAssign);
        const double mmd = transverseMMD(
            evt.genEvent.first, genClusAssign, jetClusAssign);

        // Fill out the relevant event structures
        evt.simpleDiffusionRandDist = distances.first;
        evt.simpleDiffusionJaccardDist = distances.second;
        evt.simpleDiffusionMMDDist = mmd;
        evt.simpleDiffusionClusDistReady = true;

        return true;
    }

private:
    double distanceCutoff_;
};

#endif // RANDJACCARDDMMDCALCULATOR_HH_
