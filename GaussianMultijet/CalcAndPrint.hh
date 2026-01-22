#ifndef CALCANDPRINT_HH_
#define CALCANDPRINT_HH_

#include <vector>
#include <functional>

#include "sjet/NodeVisitors.hh"

#include "AbsTextDump.hh"
#include "MemberDumper.hh"
#include "ParticleDistances.hh"
#include "matchPairToPair.hh"
#include "unpackParticle.hh"
#include "MultijetComponents.hh"

#define cap_dump_event_member(name) eventDumps_.push_back(make_MemberDumper<Event>(#name,std::mem_fn(&Event::name)))
#define cap_dump_local_member(name) localDumps_.push_back(make_MemberDumper<CalcAndPrint>(#name,std::mem_fn(&CalcAndPrint::name)))

template <class Event>
class CalcAndPrint : public frw::AbsTextDump<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsTextDump<Event> Base;

    inline CalcAndPrint(const std::string& i_label,
                        const std::string& filename)
        : Base(i_label, filename)
    {
        form_event_record();
        form_local_record();
        printLabels();
    }

    inline virtual ~CalcAndPrint() override {}

    inline virtual CalcAndPrint* clone() const override
        {return new CalcAndPrint(*this);}

    inline virtual bool analyze(const Event& evt) override
    {
        // Check that various event contents are ready
        assert(evt.isNumberValid());
        assert(evt.genEventReady);
        assert(evt.genJetsReady);
        assert(evt.diffusionSequenceReady);
        assert(evt.simpleDiffusionJetsReady);
        assert(evt.simpleDiffusionClusDistReady);

        // Match the diffusion jets with generated jets
        const Particle& refJet = evt.genJets.at(DC_REFERENCE_JET);
        const Particle& probeJet = evt.genJets.at(DC_PROBE_JET);
        unsigned j0 = 0, j1 = 0;
        nDiffusionJets = evt.simpleDiffusionJets.size();
        const ParticleDeltaR dRCalc;
        if (nDiffusionJets > 1U)
        {
            const std::pair<unsigned,unsigned>& match =
                matchPairToPair(refJet, probeJet, evt.simpleDiffusionJets, dRCalc);
            j0 = match.first;
            j1 = match.second;
        }

        // Calculate the quantities that we will be saving
        unpackParticle(refJet, refJetPt, refJetEta, refJetPhi, refJetM);
        refJetMult = evt.genEvent.second.at(DC_REFERENCE_JET);
        unpackParticle(probeJet, probeJetPt, probeJetEta, probeJetPhi, probeJetM);
        probeJetMult = evt.genEvent.second.at(DC_PROBE_JET);
        probeDrToRef = dRCalc(refJet, probeJet);

        const Particle& genDijet = refJet + probeJet;
        unpackParticle(genDijet, genDijetPt, genDijetEta, genDijetPhi, genDijetM);
        genDijetMult = refJetMult + probeJetMult;

        const Particle& j0Jet = evt.simpleDiffusionJets.at(j0);
        unpackParticle(j0Jet, j0Pt, j0Eta, j0Phi, j0M);
        {
            sjet::NodeCounter<Particle> cnt;
            evt.diffusionSequence.visitParentParticles(
                cnt, evt.simpleDiffusionNodes[j0]);
            j0Mult = cnt.result();
        }
        j0DrToRef = dRCalc(j0Jet, refJet);
        j0MsqToRef = (j0Jet - refJet).squared();

        const Particle& j1Jet = evt.simpleDiffusionJets.at(j1);
        unpackParticle(j1Jet, j1Pt, j1Eta, j1Phi, j1M);
        {
            sjet::NodeCounter<Particle> cnt;
            evt.diffusionSequence.visitParentParticles(
                cnt, evt.simpleDiffusionNodes[j1]);
            j1Mult = cnt.result();
        }
        j1DrToProbe = dRCalc(j1Jet, probeJet);
        j1MsqToProbe = (j1Jet - probeJet).squared();
        j1DrToJ0 = dRCalc(j1Jet, j0Jet);

        const Particle& diffusionDijet = j0Jet + j1Jet;
        unpackParticle(diffusionDijet, diffusionDijetPt, diffusionDijetEta,
                       diffusionDijetPhi, diffusionDijetM);
        diffusionDijetMult = j0Mult + j1Mult;

        // Calculate pileup multiplicity for diffusion jets
        diffusionPileupMult = 0;
        if (nDiffusionJets > 1U)
        {
            const unsigned jetMult = j0Mult + j1Mult;
            const unsigned nParticles = evt.diffusionSequence.nParticles();
            assert(nParticles >= jetMult);
            diffusionPileupMult = nParticles - jetMult;
        }

        // Calculate pileup scalar Pt sum for diffusion
        // jets (unclustered Pt)
        diffusionPileupScalarPtSum = 0.0;
        if (nDiffusionJets > 2U)
        {
            sjet::ScalarPtAdder<Particle> ptAdder;
            for (unsigned ijet=0; ijet<nDiffusionJets; ++ijet)
                if (ijet != j0 && ijet != j1)
                    evt.diffusionSequence.visitParentParticles(
                        ptAdder, evt.simpleDiffusionNodes[ijet]);
            diffusionPileupScalarPtSum = ptAdder.result();
        }

        // Finally, print everything
        printValues(evt);
        return true;
    }

private:
    typedef typename Event::particle_type Particle;
    typedef AbsMemberDumper<Event> AbsEventDumper;
    typedef AbsMemberDumper<CalcAndPrint> AbsLocalDumper;

    // Collections of dumpers
    std::vector<std::shared_ptr<AbsEventDumper> > eventDumps_;
    std::vector<std::shared_ptr<AbsLocalDumper> > localDumps_;

    // Various event properties to calculate and print
    //
    unsigned nDiffusionJets;  // Number of jets found. Note that if
                              // this number is less than 2, a number
                              // of variables below will not have
                              // sensible values. For making histograms,
                              // you should always select events
                              // with nDiffusionJets >= 2.

    // Characteristics of the reference jet (the one with
    // average direction (0.0, 0.0) in the eta-phi space)
    double refJetPt;          // transverse momentum
    double refJetEta;         // eta
    double refJetPhi;         // phi
    double refJetM;           // mass
    unsigned refJetMult;      // multiplicity

    // Characteristics of the probe jet
    double probeJetPt;
    double probeJetEta;
    double probeJetPhi;
    double probeJetM;
    unsigned probeJetMult;

    // Actual distance in the eta-phi space between the
    // directions of the reference and probe jets
    double probeDrToRef;

    // Characteristics of the reference and probe combined dijet system
    double genDijetPt;
    double genDijetEta;
    double genDijetPhi;
    double genDijetM;
    unsigned genDijetMult;

    // Jet 0 is matched to the reference jet
    double j0Pt;
    double j0Eta;
    double j0Phi;
    double j0M;
    unsigned j0Mult;
    double j0DrToRef;         // delta R to the reference jet 
    double j0MsqToRef;        // Squared difference between 4-momenta
                              // of this jet and reference jet

    // Jet 1 is matched to the probe jet
    double j1Pt;
    double j1Eta;
    double j1Phi;
    double j1M;
    unsigned j1Mult;
    double j1DrToProbe;       // delta R to the probe jet 
    double j1MsqToProbe;      // Squared difference between 4-momenta
                              // of this jet and probe jet
    double j1DrToJ0;          // delta R to the jet 0

    // Characteristics of the jet 0 and jet 1 combined dijet system
    double diffusionDijetPt;
    double diffusionDijetEta;
    double diffusionDijetPhi;
    double diffusionDijetM;
    unsigned diffusionDijetMult;

    unsigned diffusionPileupMult;      // Multiplicity attributed to "pileup"
                                       // particles (if more than 2 jets are
                                       // reconstructed)
    double diffusionPileupScalarPtSum; // Sum of scalar pt values for "pileup"
                                       // particles

    void printLabels() const
    {
        std::ostream& os = *Base::of_;
        os << "# ";
        printMemberDumperLabels(eventDumps_, os);
        os << ' ';
        printMemberDumperLabels(localDumps_, os);
        os << std::endl;
    }

    void printValues(const DijetEvent<Particle>& evt) const
    {
        std::ostream& os = *Base::of_;
        printMemberDumperValues(eventDumps_, evt, os);
        os << ' ';
        printMemberDumperValues(localDumps_, *this, os);
        os << std::endl;
    }

    // Decide what we are going to print from the event object
    void form_event_record()
    {
        cap_dump_event_member(number);
        cap_dump_event_member(expectedProbeEta);
        cap_dump_event_member(expectedProbePhi);
        cap_dump_event_member(simpleDiffusionRandDist);
        cap_dump_event_member(simpleDiffusionJaccardDist);
        cap_dump_event_member(simpleDiffusionMMDDist);
        cap_dump_event_member(genTotalMult);
        cap_dump_event_member(genPileupMult);
        cap_dump_event_member(genPileupScalarPtSum);
    }

    // Decide what we are going to print from this object
    void form_local_record()
    {
        cap_dump_local_member(nDiffusionJets);

        cap_dump_local_member(refJetPt);
        cap_dump_local_member(refJetEta);
        cap_dump_local_member(refJetPhi);
        cap_dump_local_member(refJetM);
        cap_dump_local_member(refJetMult);

        cap_dump_local_member(probeJetPt);
        cap_dump_local_member(probeJetEta);
        cap_dump_local_member(probeJetPhi);
        cap_dump_local_member(probeJetM);
        cap_dump_local_member(probeJetMult);
        cap_dump_local_member(probeDrToRef);

        cap_dump_local_member(genDijetPt);
        cap_dump_local_member(genDijetEta);
        cap_dump_local_member(genDijetPhi);
        cap_dump_local_member(genDijetM);
        cap_dump_local_member(genDijetMult);

        cap_dump_local_member(j0Pt);
        cap_dump_local_member(j0Eta);
        cap_dump_local_member(j0Phi);
        cap_dump_local_member(j0M);
        cap_dump_local_member(j0Mult);
        cap_dump_local_member(j0DrToRef);
        cap_dump_local_member(j0MsqToRef);

        cap_dump_local_member(j1Pt);
        cap_dump_local_member(j1Eta);
        cap_dump_local_member(j1Phi);
        cap_dump_local_member(j1M);
        cap_dump_local_member(j1Mult);
        cap_dump_local_member(j1DrToProbe);
        cap_dump_local_member(j1MsqToProbe);
        cap_dump_local_member(j1DrToJ0);

        cap_dump_local_member(diffusionDijetPt);
        cap_dump_local_member(diffusionDijetEta);
        cap_dump_local_member(diffusionDijetPhi);
        cap_dump_local_member(diffusionDijetM);
        cap_dump_local_member(diffusionDijetMult);

        cap_dump_local_member(diffusionPileupMult);
        cap_dump_local_member(diffusionPileupScalarPtSum);
    }
};

#endif // CALCANDPRINT_HH_
