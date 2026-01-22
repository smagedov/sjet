#ifndef EVENTCOUNTER_HH_
#define EVENTCOUNTER_HH_

#include <iostream>

#include "AbsFrameworkAnalyzer.hh"
#include "stringUtils.hh"

// This module can be useful if you have real filters
// in the event processing sequence or if you want to
// monitor the progress of your program by printing
// time stamps and event counts
namespace frw {
    template <class Event>
    class EventCounter : public AbsFrameworkAnalyzer<Event>
    {
    public:
        typedef AbsFrameworkAnalyzer<Event> Base;

        inline EventCounter(const std::string& i_label,
                            const unsigned printPeriod = 0U)
            : Base(i_label), count_(0), printPeriod_(printPeriod) {}
        inline virtual ~EventCounter() override {}

        inline virtual EventCounter* clone() const override
            {return new EventCounter(*this);}

        inline unsigned getCount() const {return count_;}

        inline virtual bool analyze(const Event&) override
        {
            ++count_;
            if (printPeriod_)
                if (count_ % printPeriod_ == 0U)
                {
                    if (count_ == printPeriod_)
                        std::cout << '\n';
                    std::cout << timestamp() << " : "
                              << Base::label() << " counted "
                              << count_ << " events" << std::endl;
                }
            return true;
        }

        inline virtual void endJob() const override
        {
            std::cout << Base::label() << " : counted "
                      << count_ << " events total" << std::endl;
        }

    private:
        unsigned count_;
        unsigned printPeriod_;
    };
}

#endif // EVENTCOUNTER_HH_
