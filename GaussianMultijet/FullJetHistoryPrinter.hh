#ifndef FULLJETHISTORYPRINTER_HH_
#define FULLJETHISTORYPRINTER_HH_

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
        //const int clustparts = history.size();

	of_ << node << ' ' << clus.parent1() << ' ' << clus.parent2() << ' ' << std::max(clus.dist(), 0.0) << ' ';
	printPtEtaPhiM(p, of_);
	of_ << ' ' << clus.stab();
	of_ << std::endl;

	if(clus.parent1()==-1 && clus.parent2()==-1) {
		return;
	} else {
		const cluster_type& p1 = history.at(clus.parent1());
                const cluster_type& p2 = history.at(clus.parent2());
                visit(clus.parent1(), p1.p());
                visit(clus.parent2(), p2.p());
	}

        //for (int i = 0; i < clustparts; i++) {
        //    of_ << history[i].parent1() << ' ' << history[i].parent2() << ' ' << std::max(history[i].dist(), 0.0) << ' ';
	//    printPtEtaPhiM(history[i].p(), of_);
	//    of_ << std::endl;
        //}
    }

private:
    typedef typename Event::clust_seq_type::cluster_type cluster_type;

    const Event& evt_;
    const Particle& jet_;
    std::ofstream& of_;
    bool firstCall_;
};

template <class Event>
class FullJetHistoryPrinter : public frw::AbsFrameworkAnalyzer<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkAnalyzer<Event> Base;

    inline FullJetHistoryPrinter(const std::string& i_label,
                             const std::string& outputPrefix,
                             const unsigned maxJetsToPrint,
			     const double maxDistance)
        : Base(i_label),
          prefix_(outputPrefix),
          maxJetsToPrint_(maxJetsToPrint),
	  maxDistance_(maxDistance),
          fileCount_(0)
    {
    }

    inline virtual ~FullJetHistoryPrinter() override {}

    inline virtual FullJetHistoryPrinter* clone() const override
        {return new FullJetHistoryPrinter(*this);}

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
	    of << maxDistance_ << std::endl;
            of << "node parent1 parent2 dist pt eta phi m stab"
               << std::endl;
            const Particle& thisJet = evt.simpleDiffusionJets[ijet];
            const unsigned thisNode = evt.simpleDiffusionNodes[ijet];
            JetNodePrinter prn(evt, thisJet, of);
            prn.visit(thisNode, thisJet);
            //evt.diffusionSequence.visitParentParticles(prn, thisNode);
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
    double maxDistance_;
    unsigned fileCount_;
};

#endif // FULLJETHISTORYPRINTER_HH_
