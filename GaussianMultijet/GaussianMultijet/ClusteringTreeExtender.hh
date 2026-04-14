#ifndef CLUSTERINGTREEEXTENDER_HH_
#define CLUSTERINGTREEEXTENDER_HH_

#include "AbsFrameworkModule.hh"
#include "sjet/DistanceCalculator.hh"

template <class Event>
class ClusteringTreeExtender : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;
    typedef std::vector<typename Event::particle_type> ParticleVector;

    inline ClusteringTreeExtender(const std::string& i_label)
        : Base(i_label)
    {
    }

    inline virtual ~ClusteringTreeExtender() override {}

    inline virtual ClusteringTreeExtender* clone() const override
        {return new ClusteringTreeExtender(*this);}

    inline virtual bool process(Event& evt) override
    {
        assert(evt.diffusionSequenceReady);
        assert(!evt.copySequenceReady);

        ParticleVector extendedInfo;
        extendedInfo.reserve(evt.diffusionSequence.nClusters());

        evt.copySequence = event_type::extended_seq_type(
            evt.diffusionSequence, sjet::DummyCalculator(), extendedInfo);
        evt.copySequenceReady = true;
        return true;
    }

private:
};

#endif // CLUSTERINGTREEEXTENDER_HH_
