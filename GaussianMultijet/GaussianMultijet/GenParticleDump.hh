#ifndef GENPARTICLEDUMP_HH_
#define GENPARTICLEDUMP_HH_

#include <cassert>
#include <set>

#include "AbsTextDump.hh"
#include "printPtEtaPhi.hh"

template <class Event>
class GenParticleDump : public frw::AbsTextDump<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsTextDump<Event> Base;

    // If this constructor is used, all events will be dumped
    inline GenParticleDump(const std::string& i_label,
                           const std::string& filename)
        : Base(i_label, filename),
          dumpAll_(true),
          firstDump_(true)
    {
    }

    // Dump one particular event
    inline GenParticleDump(const std::string& i_label,
                           const std::string& filename,
                           const unsigned long eventToDump)
        : Base(i_label, filename),
          dumpAll_(false),
          firstDump_(true)
    {
        eventSet_.insert(eventToDump);
    }

    // Dump events specified in the input collection
    // (normally, vector or set)
    template<class Collection>
    inline GenParticleDump(const std::string& i_label,
                           const std::string& filename,
                           const Collection& eventsToDump)
        : Base(i_label, filename),
          eventSet_(eventsToDump.begin(), eventsToDump.end()),
          dumpAll_(false),
          firstDump_(true)
    {
    }

    inline virtual ~GenParticleDump() override {}

    inline virtual GenParticleDump* clone() const override
        {return new GenParticleDump(*this);}

    inline virtual bool analyze(const Event& evt) override
    {
        assert(evt.isNumberValid());
        assert(evt.genEventReady);

        // Should we dump this event?
        if (shouldDumpThis(evt.number()))
        {
            if (firstDump_)
            {
                *Base::of_ << "# pt eta phi" << std::endl;
                firstDump_ = false;
            }

            *Base::of_ << "# Event " << evt.number() << std::endl;
            for (const Particle& p : evt.genEvent.first)
            {
                printPtEtaPhi(p, *Base::of_);
                *Base::of_ << std::endl;
            }
        }

        return true;
    }

private:
    typedef typename Event::particle_type Particle;

    inline bool shouldDumpThis(const unsigned long evtNumber) const
    {
        if (dumpAll_)
            return true;
        else
            return !(eventSet_.find(evtNumber) == eventSet_.end());
    }

    std::set<unsigned long> eventSet_;
    bool dumpAll_;
    bool firstDump_;
};

#endif // GENPARTICLEDUMP_HH_
