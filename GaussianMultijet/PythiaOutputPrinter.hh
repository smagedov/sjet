#ifndef PYTHIAOUTPUTPRINTER_HH_
#define PYTHIAOUTPUTPRINTER_HH

#include <cassert>
#include <set>

#include "AbsTextDump.hh"
#include "printPtEtaPhi.hh"
#include "ParticleDistances.hh"

template <class Event>
class PythiaOutputPrinter : public frw::AbsTextDump<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsTextDump<Event> Base;

    // If this constructor is used, all events will be dumped
    inline PythiaOutputPrinter(const std::string& i_label,
                           const std::string& outname)
        : Base(i_label, outname),
          firstDump_(true)
    {
    }

    inline virtual ~PythiaOutputPrinter() override {}

    inline virtual PythiaOutputPrinter* clone() const override
        {return new PythiaOutputPrinter(*this);}

    inline virtual bool analyze(const Event& evt) override
    {
        assert(evt.genEventReady);
	std::vector<Pythia8::Particle> finalParticles;
	finalParticles.reserve(evt.pythiaEvent->size());
	
	// Loop through particles in the event
	for (int i = 0; i < evt.pythiaEvent->size(); ++i) {
		const Pythia8::Particle& p = (*evt.pythiaEvent)[i];
		if (p.isFinal()) {
			finalParticles.push_back(p);
		}
	}
	
	int Jets = 0;
	for (long unsigned int i=0; i<evt.genClusters.size(); ++i) {
		if (!evt.genClusters[i].empty()) {
			Jets = Jets + 1;
			int partCount = evt.genClusters[i].size();
			*Base::of_ << "Jet #" << Jets << " Particle Count: " << partCount << std::endl;
		}
	}
	
	*Base::of_ << "nJets: " << Jets << std::endl;;
	*Base::of_ << "id" << " " << "pT" << " " << "eta" << " " << "phi" << " " << "m" << std::endl;
	
	for (Pythia8::Particle p : finalParticles) {
		*Base::of_ << p.id()     << " "
		           << p.pT()     << " "
		           << p.eta()     << " "
		           << p.phi()     << " "
		           << p.m()      << std::endl;
	}
	
	return true;
    }

private:
    typedef typename Event::particle_type Particle;

    bool firstDump_;
};

#endif // PYTHIAOUTPUTPRINTER_HH_
