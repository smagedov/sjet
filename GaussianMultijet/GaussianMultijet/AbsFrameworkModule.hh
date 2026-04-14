#ifndef FRW_ABSFRAMEWORKMODULE_HH_
#define FRW_ABSFRAMEWORKMODULE_HH_

#include <string>

namespace frw {
    template <class Event>
    class AbsFrameworkModule
    {
    public:
        typedef Event event_type;

        // Every module must have a label (for identification
        // in various printouts, etc). All module labels in
        // a module sequence must be unique.
        inline AbsFrameworkModule(const std::string& i_label)
            : label_(i_label) {}

        inline virtual ~AbsFrameworkModule() {}

        // "Virtual copy constructor"
        virtual AbsFrameworkModule* clone() const = 0;

        // Return the module label
        inline const std::string& label() const {return label_;}

        // The following function must return "true" if processing
        // of this event by the module sequence should continue past
        // this module, "false" otherwise
        virtual bool process(Event& evt) = 0;

        // The following function will be called once at the end of
        // event processing. Override this method to save objects,
        // print summaries, etc. Note that the original module is
        // typically destroyed when a copy is placed into the sequence
        // of modules, so don't save/print things in the destructor as
        // this will usually happen twice.
        inline virtual void endJob() const {}

    private:
        std::string label_;
    };
}

#endif // FRW_ABSFRAMEWORKMODULE_HH_
