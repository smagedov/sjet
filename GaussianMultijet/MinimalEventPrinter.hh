#ifndef MINIMALEVENTPRINTER_HH_
#define MINIMALEVENTPRINTER_HH_

#include <cassert>
#include <limits>
#include <set>

#include "AbsTextDump.hh"
#include "printPtEtaPhi.hh"

template <class Event>
struct MinimalEventPrinter : public frw::AbsTextDump<Event>
{
    typedef Event event_type;
    typedef frw::AbsTextDump<Event> Base;

    // If this constructor is called without specifying the event
    // number to dump, all events will be dumped
    inline MinimalEventPrinter(const std::string& i_label,
                               const std::string& filename)
        : Base(i_label, filename) {}

    inline virtual ~MinimalEventPrinter() override {}

    inline virtual MinimalEventPrinter* clone() const override
        {return new MinimalEventPrinter(*this);}

    inline virtual bool analyze(const Event& evt) override
    {
        // Make sure that the portion of the data we are
        // going to work with is ready
        assert(evt.randomDatumReady);

        // Write this data into a file
        *Base::of_ << "Event " << evt.number()
                   << " : " << evt.randomDatum << std::endl;

        return true;
    }
};

#endif // MINIMALEVENTPRINTER_HH_
