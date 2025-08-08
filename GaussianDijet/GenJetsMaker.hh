#ifndef GENJETSMAKER_HH_
#define GENJETSMAKER_HH_

#include "AbsFrameworkModule.hh"

template <class Event>
struct GenJetsMaker : public frw::AbsFrameworkModule<Event>
{
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline GenJetsMaker(const std::string& i_label)
        : Base(i_label) {}

    inline virtual ~GenJetsMaker() override {}

    inline virtual GenJetsMaker* clone() const override
        {return new GenJetsMaker(*this);}

    inline virtual bool process(Event& evt) override
    {
        // Make sure that we have the input data structures ready
        assert(evt.genEventReady);

        // Make sure that we are not overwriting a part of the
        // event already made by some other module
        assert(!evt.genJetsReady);

        // Sum up the particles which make generated jets
        evt.genJets.clear();
        evt.genJets.reserve(evt.genEvent.second.size());
        unsigned imin = 0;
        for (unsigned blockLength : evt.genEvent.second)
        {
            Particle jet;
            const unsigned imax = imin + blockLength;
            for (unsigned i=imin; i<imax; ++i)
                jet += evt.genEvent.first.at(i);
            imin = imax;
            evt.genJets.push_back(jet);
        }

        evt.genJetsReady = true;
        return true;
    }

private:
    typedef typename Event::particle_type Particle;
};

#endif // GENJETSMAKER_HH_
