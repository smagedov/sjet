#ifndef FRW_ABSFRAMEWORKANALYZER_HH_
#define FRW_ABSFRAMEWORKANALYZER_HH_

#include "AbsFrameworkModule.hh"

// AbsFrameworkAnalyzer is accessing event by const reference
// in its "analyze" function, other than that things are
// similar to AbsFrameworkModule.
namespace frw {
    template <class Event>
    class AbsFrameworkAnalyzer : public AbsFrameworkModule<Event>
    {
    public:
        typedef Event event_type;
        typedef AbsFrameworkModule<Event> Base;

        inline AbsFrameworkAnalyzer(const std::string& i_label)
            : Base(i_label) {}

        inline virtual ~AbsFrameworkAnalyzer() override {}

        // "Virtual copy constructor"
        virtual AbsFrameworkAnalyzer* clone() const = 0;

        // The following function must return "true" if processing
        // of this event by the module sequence should continue past
        // this module, "false" otherwise
        virtual bool analyze(const Event& evt) = 0;

        // The following function will be called once at the end of
        // event processing. Override this method to save objects,
        // print summaries, etc. Note that the original module is
        // typically destroyed when a copy is placed into the sequence
        // of modules, so don't save/print things in the destructor as
        // this will usually happen twice.
        inline virtual void endJob() const override {}

    private:
        inline bool process(Event& evt) override
            {return analyze(evt);}
    };
}

#endif // FRW_ABSFRAMEWORKANALYZER_HH_
