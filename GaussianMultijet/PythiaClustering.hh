#ifndef PYTHIACLUSTERING_HH_
#define PYTHIACLUSTERING_HH_

#include <cassert>
#include <cfloat>
#include <set>

#include "AbsTextDump.hh"
#include "printPtEtaPhi.hh"
#include "ParticleDistances.hh"
#include "permutation.hh"

template <class Event>
class PythiaClustering : public frw::AbsTextDump<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsTextDump<Event> Base;

    // If this constructor is used, all events will be dumped
    inline PythiaClustering(const std::string& i_label,
                           const std::string& filename)
        : Base(i_label, filename),
	  firstDump_(true)
    {
    }

    inline virtual ~PythiaClustering() override {}

    inline virtual PythiaClustering* clone() const override
        {return new PythiaClustering(*this);}

    inline virtual bool analyze(const Event& evt) override
    {
        assert(evt.genEventReady);
	assert(evt.simpleDiffusionJetsReady);

	const std::vector<unsigned>& genClusAssign =
        genClusterAssignments(evt.genEvent.second);        

        std::vector<unsigned> jetClusAssign =
        evt.diffusionSequence.clusterAssignments(evt.simpleDiffusionJetsDistCutoff);
	
	// Find which node in the clustering tree corresponds
	// closest to the reference jet
	std::vector<unsigned> jetComponents;
	unsigned comp = 0;
	for (unsigned i=0; i<evt.genClusters.size(); ++i) {
                if (!evt.genClusters[i].empty()) {
			jetComponents.push_back(comp);
			comp = comp + 1;
		}
	}

	const long unsigned nDiffusionJets = evt.simpleDiffusionJets.size();
        const long unsigned nGenJets = evt.genJets.size();

	std::cout << "GenJets: " << nGenJets << " ClusteredJets: " << nDiffusionJets << std::endl;

        auto perms = generatePermutations(nDiffusionJets, nGenJets);

	std::cout << "Total Permutations: " << perms.size() << std::endl;

	const ParticleDeltaR dRcalculator;
	double minDist = DBL_MAX;
	std::vector<int> bestAssign(nGenJets);
	for (unsigned i=0; i<perms.size(); ++i) {
		double totDist = 0;
		for (unsigned j=0; j<perms[i].size(); ++j) {
			const Particle& genJet = evt.genJets.at(j);
			const Particle& difJet = evt.simpleDiffusionJets.at(perms[i][j]);
			double tmp=dRcalculator(genJet, difJet);
			totDist = totDist + tmp;
		}
		if (totDist < minDist)
		{
			bestAssign = perms[i];
			minDist = totDist;
		}
	}

	for (unsigned i=0; i<bestAssign.size(); ++i) {
		std::cout << bestAssign[i] << " ";
	}
	std::cout << std::endl;

        const unsigned nParticles = jetClusAssign.size();

	std::cout << "Simple Diffusion Nodes Size: " << evt.simpleDiffusionNodes.size() << " nParticles: " << nParticles  << std::endl;

        std::vector<unsigned> clustNumberAssign;
	//std::vector<unsigned> genClusAssign;
	for (unsigned i=0; i<nParticles; ++i) {
		clustNumberAssign.push_back(bestAssign[genClusAssign[i]]);
	}

        if (firstDump_)
        {
            *Base::of_ << "# pt eta phi true_cluster_assignment clustering_cluster_assignment" << std::endl;
	    firstDump_ = false;
        }

        *Base::of_ << "# Event " << evt.number() << std::endl;

	for (unsigned i=0; i<nParticles; ++i)
	{
            printPtEtaPhi(evt.genEvent.first[i], *Base::of_);
	    *Base::of_ << ' ' << genClusAssign[i] << ' ' << clustNumberAssign[i];
            *Base::of_ << std::endl;
        }

        return true;
    }

private:
    typedef typename Event::particle_type Particle;

    bool firstDump_;
};

#endif // PYTHIACLUSTERING_HH_
