#ifndef MINIMALEVENTMAKER_HH_
#define MINIMALEVENTMAKER_HH_

#include <random>
#include <cassert>

#include "AbsFrameworkModule.hh"

template <class Event, class Rng>
class MinimalEventMaker : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef Rng rng_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline MinimalEventMaker(const std::string& i_label,
                             Rng& gen, const double maxValue)
        : Base(i_label), gen_(gen), uni_(0.0, maxValue) {}

    inline virtual ~MinimalEventMaker() override {}

    inline virtual MinimalEventMaker* clone() const override
        {return new MinimalEventMaker(*this);}

    inline virtual bool process(Event& evt) override
    {
        // Make sure that the event has been initialized
        assert(evt.isNumberValid());

        // Make sure that we are not overwriting existing data
        assert(!evt.randomDatumReady);

        // Fill out this module's portion of event data
        evt.randomDatum = uni_(gen_);

        // Note that this portion of event is ready
        evt.randomDatumReady = true;

        // Return allowing other modules to proceed
        return true;
    }

private:
    Rng& gen_;
    std::uniform_real_distribution<> uni_;
};

#endif // MINIMALEVENTMAKER_HH_
