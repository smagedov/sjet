#ifndef PYTHIACLUSTERING_HH_
#define PYTHIACLUSTERING_HH_

#include <cassert>
#include <set>

#include "AbsTextDump.hh"
#include "printPtEtaPhi.hh"
#include "ParticleDistances.hh"

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

	std::vector<unsigned> clusComponents;
	double totDist = 0;
	const unsigned nDiffusionJets = evt.simpleDiffusionJets.size();
        const ParticleDeltaR dRcalculator;

	for (unsigned i=0; i<jetComponents.size(); ++i) {
		const Particle& refJet = evt.genJets.at(jetComponents[i]);
		unsigned closesti=0;
		double minDist=dRcalculator(refJet, evt.simpleDiffusionJets[closesti]);
		for (unsigned j=1; j<nDiffusionJets; j++)
		{
			const double tmp = dRcalculator(refJet, evt.simpleDiffusionJets[j]);
			if (tmp < minDist)
			{
				closesti = j;
				minDist = tmp;
			}
		}
		clusComponents.push_back(closesti);
		totDist = totDist + minDist;
	}

	std::vector<unsigned> closestNode;
	for (unsigned i=0; i<clusComponents.size(); ++i) {
		closestNode.push_back(evt.simpleDiffusionNodes[clusComponents[i]]);
	}

	// Only use the diffusion clustering assignments to the
        // two leading pt jets. Assume everything else is pileup
        // or unclustered.
	const unsigned nParticles = jetClusAssign.size();

        if (nDiffusionJets > jetComponents.size())
	{
            const unsigned unusedNode = evt.diffusionSequence.clustHist().size();

            for (unsigned i=0; i<nParticles; ++i)
            {
                const unsigned a = jetClusAssign[i];
                if (std::find(jetComponents.begin(), jetComponents.end(), a) == jetComponents.end())
                    jetClusAssign[i] = unusedNode;
            }
        }

	std::vector<unsigned> clustNumberAssign(nParticles);

	for (unsigned i=0; i<nParticles; ++i)
	{
		for (unsigned j=0; j<closestNode.size(); ++j) {
			std::cout << "Pong" << std::endl;
			if (jetClusAssign[i] == closestNode[j]) {
				std::cout << "Ping" << std::endl;
				clustNumberAssign[i] = jetComponents[j];
			}
		}
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
