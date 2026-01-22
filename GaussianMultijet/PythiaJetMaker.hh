#ifndef PYTHIAJETMAKER_HH_
#define PYTHIAJETMAKER_HH_

#include <cassert>
#include <algorithm>
#include <utility>
#include <vector>

#include "AbsFrameworkModule.hh"

template <class Event>
class PythiaJetMaker : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline PythiaJetMaker(const std::string& i_label,
                          const double distanceCutoff)
        : Base(i_label), distanceCutoff_(distanceCutoff)
    {
        assert(distanceCutoff_ > 0.0);
    }

    inline virtual ~PythiaJetMaker() override {}

    inline virtual PythiaJetMaker* clone() const override
        {return new PythiaJetMaker(*this);}

    inline virtual bool process(Event& evt) override
    {
        // Make sure that we have the input data structures ready
        assert(evt.diffusionSequenceReady);

        // Make sure that we are not overwriting a part of the
        // event already made by some other module
        assert(!evt.simpleDiffusionJetsReady);

        // Extract the cluster indices from the clustering history
        const std::vector<unsigned>& indices =
            evt.diffusionSequence.clusterIndices(distanceCutoff_);

        // Sort these clusters in the order of decreasing pt
        const unsigned nClus = indices.size();
        std::vector<std::pair<double, unsigned> > ptIndices;
        ptIndices.reserve(nClus);
        const auto& clusHistory = evt.diffusionSequence.clustHist();
        for (unsigned i=0; i<nClus; ++i)
            ptIndices.emplace_back(clusHistory[indices[i]].p().pt(), indices[i]);
        std::sort(ptIndices.begin(), ptIndices.end(),
                  std::greater<std::pair<double, unsigned> >());

        // Now, fill the relevant event structures
        evt.simpleDiffusionJets.clear();
        evt.simpleDiffusionJets.reserve(nClus);
        evt.simpleDiffusionNodes.clear();
        evt.simpleDiffusionNodes.reserve(nClus);
        for (unsigned i=0; i<nClus; ++i)
        {
            const unsigned node = ptIndices[i].second;
            evt.simpleDiffusionNodes.push_back(node);
            evt.simpleDiffusionJets.push_back(clusHistory[node].p());
        }
        evt.simpleDiffusionJetsDistCutoff = distanceCutoff_;
        evt.simpleDiffusionJetsReady = true;

        return true;
    }

private:
    double distanceCutoff_;
};

#endif // PYTHIAJETMAKER_HH_
