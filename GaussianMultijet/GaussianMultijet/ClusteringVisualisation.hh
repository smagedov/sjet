#ifndef CLUSTERINGVISUALISATION_HH_
#define CLUSTERINGVISUALISATION_HH_

#include <cassert>
#include <set>

#include "AbsTextDump.hh"
#include "printPtEtaPhi.hh"
#include "ParticleDistances.hh"

template <class Event>
class ClusteringVisualisation : public frw::AbsTextDump<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsTextDump<Event> Base;

    // If this constructor is used, all events will be dumped
    inline ClusteringVisualisation(const std::string& i_label,
                           const std::string& filename)
        : Base(i_label, filename),
	  firstDump_(true)
    {
    }

    inline virtual ~ClusteringVisualisation() override {}

    inline virtual ClusteringVisualisation* clone() const override
        {return new ClusteringVisualisation(*this);}

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
	const Particle& refJet = evt.genJets.at(DC_REFERENCE_JET);
	const unsigned nDiffusionJets = evt.simpleDiffusionJets.size();
	const ParticleDeltaR dRcalculator;
	unsigned closestRefi=0;
	double minDist=dRcalculator(refJet, evt.simpleDiffusionJets[closestRefi]);
	for (unsigned i=1; i<nDiffusionJets; i++)
	{
		const double tmp = dRcalculator(refJet, evt.simpleDiffusionJets[i]);
		if (tmp < minDist)
		{
			closestRefi = i;
			minDist = tmp;
		}
	}

	const unsigned closestRefNode = evt.simpleDiffusionNodes[closestRefi];
	unsigned closestProbi=0;
	if (closestRefi == 0) {
		if (nDiffusionJets > 1) {
			closestProbi = 1;
		} else {
		// What to do if only 1 Jet is reconstructed
		//
		}
	}
	const unsigned closestProbNode = evt.simpleDiffusionNodes[closestProbi];

	// Only use the diffusion clustering assignments to the
        // two leading pt jets. Assume everything else is pileup
        // or unclustered.
	const unsigned nParticles = jetClusAssign.size();

        if (nDiffusionJets > DC_N_EVENT_COMPONENTS)
        {
            const unsigned node0 = evt.simpleDiffusionNodes[0];
            const unsigned node1 = evt.simpleDiffusionNodes[1];
            const unsigned unusedNode = evt.diffusionSequence.clustHist().size();

            for (unsigned i=0; i<nParticles; ++i)
            {
                const unsigned a = jetClusAssign[i];
                if (a != node0 && a != node1)
                    jetClusAssign[i] = unusedNode;
            }
        }

	std::vector<unsigned> clustNumberAssign(nParticles);

	for (unsigned i=0; i<nParticles; ++i)
	{
		if (jetClusAssign[i] == closestRefNode) {
			clustNumberAssign[i] = DC_REFERENCE_JET;
		} else if (jetClusAssign[i] == closestProbNode) {
			clustNumberAssign[i] = DC_PROBE_JET;
		} else {
			clustNumberAssign[i] = DC_PILEUP;
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

#endif // CLUSTERINGVISUALISATION_HH_
