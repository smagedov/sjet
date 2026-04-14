#ifndef FRW_FRAMEWORKANALYSISSEQUENCE_HH_
#define FRW_FRAMEWORKANALYSISSEQUENCE_HH_

#include <memory>
#include <vector>
#include <stdexcept>
#include <sstream>

#include "AbsFrameworkModule.hh"

namespace frw {
    template <class Event>
    class FrameworkAnalysisSequence
    {
    public:
        inline FrameworkAnalysisSequence() {}

        inline unsigned nModules() const {return modules_.size();}

        inline void add(const AbsFrameworkModule<Event>& i_module)
        {
            if (findByLabel(i_module.label()))
            {
                std::ostringstream os;
                os << "In frw::FrameworkAnalysisSequence::add : "
                   << "duplicate module label \""
                   << i_module.label() << "\" is not alowed";
                throw std::invalid_argument(os.str());
            }
            modules_.emplace_back(i_module.clone());
        }

        inline AbsFrameworkModule<Event>& module(const unsigned i)
            {return *modules_.at(i);}

        inline const AbsFrameworkModule<Event>& module(const unsigned i) const
            {return *modules_.at(i);}

        inline std::vector<std::string> labels() const
        {
            std::vector<std::string> result;
            const unsigned sz = modules_.size();
            result.reserve(sz);
            for (unsigned imod=0; imod<sz; ++imod)
                result.push_back(modules_[imod]->label());
            return result;
        }

        // The following function returns nullptr in case a module
        // with this label is not found
        inline const AbsFrameworkModule<Event>* findByLabel(
            const std::string& l) const
        {
            const unsigned sz = modules_.size();
            for (unsigned imod=0; imod<sz; ++imod)
                if (l == modules_[imod]->label())
                    return &*modules_[imod];
            return nullptr;
        }

        inline void run(Event& evt)
        {
            const unsigned sz = modules_.size();
            for (unsigned imod=0; imod<sz; ++imod)
                if (!modules_[imod]->process(evt))
                    break;
        }

        inline void endJob() const
        {
            const unsigned sz = modules_.size();
            for (unsigned imod=0; imod<sz; ++imod)
                modules_[imod]->endJob();
        }

    private:
        // Objects of this class should not be copied
        FrameworkAnalysisSequence(const FrameworkAnalysisSequence&) = delete;
        FrameworkAnalysisSequence& operator=(const FrameworkAnalysisSequence&) = delete;

        std::vector<std::shared_ptr<AbsFrameworkModule<Event> > > modules_;
    };
}

#endif // FRW_FRAMEWORKANALYSISSEQUENCE_HH_
