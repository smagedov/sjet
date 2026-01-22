#ifndef FRW_ABSTEXTDUMP_HH_
#define FRW_ABSTEXTDUMP_HH_

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "AbsFrameworkAnalyzer.hh"

// AbsTextDump is an analyser which writes its output into some
// text file
namespace frw {
    template <class Event>
    class AbsTextDump : public AbsFrameworkAnalyzer<Event>
    {
    public:
        typedef Event event_type;
        typedef AbsFrameworkAnalyzer<Event> Base;

        inline AbsTextDump(const std::string& i_label,
                           const std::string& filename,
                           const unsigned precision = 12U)
            : Base(i_label),
              outfile_(filename),
              of_(new std::ofstream(outfile_))
        {
            if (of_->is_open())
                of_->precision(precision);
            else
            {
                std::ostringstream os;
                os << "In frw::AbsTextDump constructor : "
                   << "failed to open file \""
                   << outfile_ << '"';
                throw std::invalid_argument(os.str());
            }
        }

        inline virtual ~AbsTextDump() override {}

        virtual AbsTextDump* clone() const = 0;
        virtual bool analyze(const Event& evt) = 0;

        inline virtual void endJob() const override
        {
            of_->close();
            std::cout << Base::label() << " : wrote file \""
                      << outfile_ << '"' << std::endl;
        }

    protected:
        // std::ofstream does not have a copy constructor,
        // so we are using a shared_ptr here in order to
        // enable default copy constructor generation
        std::string outfile_;
        std::shared_ptr<std::ofstream> of_;
    };
}

#endif // FRW_ABSTEXTDUMP_HH_
