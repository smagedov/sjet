#ifndef RUNDIFFUSIONCLUSTERING_HH_
#define RUNDIFFUSIONCLUSTERING_HH_

#include "AbsFrameworkModule.hh"

template <class Event>
class RunDiffusionClustering : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline RunDiffusionClustering(const std::string& i_label,
                                  const double maxDistance)
        : Base(i_label), maxDist_(maxDistance) {}

    inline virtual ~RunDiffusionClustering() override {}

    inline virtual RunDiffusionClustering* clone() const override
        {return new RunDiffusionClustering(*this);}

    inline virtual bool process(Event& evt) override
    {
        // Make sure that we have the input data structures ready
        assert(evt.genEventReady);

        // Make sure that we are not overwriting a part of the
        // event already made by some other module
        assert(!evt.diffusionSequenceReady);

        // Sum up the particles which make generated jets
        evt.diffusionSequence.init(evt.genEvent.first);
        evt.diffusionSequence.run(maxDist_);
        evt.diffusionSequenceReady = true;

        return true;
    }

private:
    double maxDist_;
};

#endif // RUNDIFFUSIONCLUSTERING_HH_
