#ifndef MULTIJETMAKER_HH_
#define MULTIJETMAKER_HH_

#include <cmath>
#include <random>

#include "AbsFrameworkModule.hh"
#include "MultiParticleCollectionMaker.hh"
#include "GaussianJetV1.hh"
#include "ExponentialPileup.hh"
#include "filterByPt.hh"

template <class Event, class Rng>
class MultijetMaker : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef Rng rng_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline MultijetMaker(const std::string& i_label,
                      Rng& gen,
                      const double refJetWidth,
                      const double refJetAvePt,
                      const double refJetMult,
                      const double probeJetWidth,
                      const double probeJetAvePt,
                      const double probeJetMult,
                      const double pileupEtaRange,
                      const double pileupAvePt,
                      const double puleupMult,
                      const double deltaR,
                      const double ptMin,
                      const bool randomizeProbeJetDirection,
                      const bool randomizeMultiplicity,
		      const unsigned nJets)
        : Base(i_label), gen_(gen), uniCircle_(-M_PI, M_PI), deltaR_(deltaR),
          ptMin_(ptMin), randomizeProbe_(randomizeProbeJetDirection)
    {
        const MyJet refJet(0.0, 0.0, refJetWidth, refJetAvePt,
                           refJetMult, randomizeMultiplicity);
        eventBuilder_.add(refJet);

	// Define eta and phi ranges
        // double eta_min = 0;
        // double eta_max = 10;

        // double phi_min = -2.5;
        // double phi_max = 2.5;

	// Set up random number generators
        // std::uniform_real_distribution<double> eta_dist(eta_min, eta_max);
        // std::uniform_real_distribution<double> phi_dist(phi_min, phi_max);

	for(unsigned i = 0; i < nJets; i ++) {
		const double etaPhiAngle = randomizeProbe_ ? uniCircle_(gen_) : static_cast<double>(i);
                double expectedProbeEta = deltaR_*cos(etaPhiAngle);
                double expectedProbePhi = deltaR_*sin(etaPhiAngle);

		const MyJet probJet(expectedProbeEta, expectedProbePhi, probeJetWidth, probeJetAvePt,
				probeJetMult, randomizeMultiplicity);
		eventBuilder_.add(probJet);
	}

        const MyPileup pileup(pileupEtaRange, pileupAvePt,
                              puleupMult, randomizeMultiplicity);
        eventBuilder_.add(pileup);

        assert(eventBuilder_.nMakers() == DC_N_EVENT_COMPONENTS);
    }

    inline virtual ~MultijetMaker() override {}

    inline virtual MultijetMaker* clone() const override
        {return new MultijetMaker(*this);}

    inline virtual bool process(Event& evt) override
    {
        // Make sure that the event has been initialized
        assert(evt.isNumberValid());

        // Make sure that we are not overwriting existing data
        assert(!evt.genEventReady);

        // If requested, randomize the location of the probe jet.
        // Its nominal center will remain at the given deltaR
        // from the reference jet.
        // const double etaPhiAngle = randomizeProbe_ ? uniCircle_(gen_) : 0.0;
        // evt.expectedProbeEta = deltaR_*cos(etaPhiAngle);
        // evt.expectedProbePhi = deltaR_*sin(etaPhiAngle);
	// eventBuilder_.maker(DC_PROBE_JET).setLocation(
	//		evt.expectedProbeEta, evt.expectedProbePhi);

        // Generate the initial event (not yet filtered by Pt)
        const GenEvent& genEvent0 = eventBuilder_.make(gen_);

        // Filter particles by their Pt
        filterByPt(genEvent0, &evt.genEvent, ptMin_);

        evt.genEventReady = true;
        return true;
    }

private:
    typedef typename Event::particle_type MyParticle;
    typedef GaussianJetV1<MyParticle,Rng> MyJet;
    typedef ExponentialPileup<MyParticle,Rng> MyPileup;
    typedef MultiParticleCollectionMaker<MyParticle,Rng> MyEventBuilder;
    typedef std::pair<std::vector<MyParticle>, std::vector<unsigned> > GenEvent;

    Rng& gen_;
    MyEventBuilder eventBuilder_;
    std::uniform_real_distribution<> uniCircle_;
    double deltaR_;
    double ptMin_;
    bool randomizeProbe_;
};

#endif // MULTIJETMAKER_HH_
