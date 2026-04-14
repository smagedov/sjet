#ifndef JETHISTORYPRINTER_HH_
#define JETHISTORYPRINTER_HH_

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "sjet/AbsNodeVisitor.hh"

#include "AbsFrameworkAnalyzer.hh"
#include "stringUtils.hh"
#include "printPtEtaPhi.hh"
#include "ParticleComparators.hh"
#include "mathUtils.hh"

template <class Event>
class JetNodePrinter : public sjet::AbsNodeVisitor<typename Event::particle_type>
{
public:
    typedef typename Event::particle_type Particle;

    inline JetNodePrinter(const Event& evt, const Particle& jet,
                          std::ofstream& of)
        : evt_(evt), jet_(jet), of_(of), firstCall_(true) {}

    inline virtual ~JetNodePrinter() {}

    virtual void visit(const unsigned node, const Particle& p)
    {
        const std::vector<cluster_type>& history = evt_.diffusionSequence.clustHist();
        const cluster_type& clus = history.at(node);
        const double dist = std::max(clus.dist(), 0.0);
        const double sPt = clus.scalarPtSum();

        of_ << dist << ' ' << std::max(clus.maxDistanceSoFar(), 0.0) << ' ' << sPt << ' ';
        printPtEtaPhiM(p, of_);

        const double pt = p.pt();
        for (int power=-2; power<0; ++power)
        {
            const double tmp = safepow(dist, power);
            of_ << ' ' << pt*tmp << ' ' << sPt*tmp;
        }

        // Pt fraction of this particle in the daughter and in
        // the overall jet, delta R to dau and to the overall jet
        if (firstCall_)
        {
            of_ << " 1.0 1.0 0.0 0.0";
            firstCall_ = false;
        }
        else
        {
            const cluster_type& dau = history.at(clus.daughter());
            of_ << ' ' << pt/dau.p().pt() << ' ' << pt/jet_.pt();
            of_ << ' ' << geom3::deltaR(p.momentum(), dau.p().momentum());
            of_ << ' ' << geom3::deltaR(p.momentum(), jet_.momentum());
        }

        of_ << std::endl;
    }

private:
    typedef typename Event::clust_seq_type::cluster_type cluster_type;

    const Event& evt_;
    const Particle& jet_;
    std::ofstream& of_;
    bool firstCall_;
};

template <class Event>
class JetHistoryPrinter : public frw::AbsFrameworkAnalyzer<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkAnalyzer<Event> Base;

    inline JetHistoryPrinter(const std::string& i_label,
                             const std::string& outputPrefix,
                             const unsigned maxJetsToPrint)
        : Base(i_label),
          prefix_(outputPrefix),
          maxJetsToPrint_(maxJetsToPrint),
          fileCount_(0)
    {
    }

    inline virtual ~JetHistoryPrinter() override {}

    inline virtual JetHistoryPrinter* clone() const override
        {return new JetHistoryPrinter(*this);}

    inline virtual bool analyze(const Event& evt) override
    {
        assert(evt.diffusionSequenceReady);
        assert(evt.simpleDiffusionJetsReady);

        const unsigned nJetNodes = evt.simpleDiffusionNodes.size();
        assert(nJetNodes == evt.simpleDiffusionJets.size());
        const unsigned nPrint = std::min(maxJetsToPrint_, nJetNodes);
        for (unsigned ijet=0; ijet<nPrint; ++ijet)
        {
            const std::string& fname = concatStringsAndInts(
                prefix_, evt.number(), ijet, ".txt");
            std::ofstream of(fname);
            if (!of.is_open())
            {
                std::ostringstream os;
                os << "In JetHistoryPrinter::analyze : failed to ope file \""
                   << fname << '"';
                throw std::invalid_argument(os.str());
            }
            of.precision(12);
            of << "# parent1 parent2 dist maxDistSoFar sPt pt eta phi m"
               << " ptdistm2 sPtdistm2 ptdistm1 sPtistdm1"
               << " ptFracDau ptFracJet dRDau dRJet"
               << std::endl;
            const Particle& thisJet = evt.simpleDiffusionJets[ijet];
            const unsigned thisNode = evt.simpleDiffusionNodes[ijet];
            JetNodePrinter prn(evt, thisJet, of);
            prn.visit(thisNode, thisJet);
            evt.diffusionSequence.visitLargerAncestors(prn, thisNode, LessByPt<Particle>());
            ++fileCount_;
        }

        return true;
    }

    inline virtual void endJob() const override
    {
        std::cout << Base::label() << " : wrote "
                  << fileCount_ << " files with prefix \""
                  << prefix_ << '"' << std::endl;
    }

private:
    typedef typename Event::particle_type Particle;

    std::string prefix_;
    unsigned maxJetsToPrint_;
    unsigned fileCount_;
};

#endif // JETHISTORYPRINTER_HH_
