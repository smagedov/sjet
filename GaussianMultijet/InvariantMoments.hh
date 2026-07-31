#ifndef INVARIANTMOMENTS_HH_
#define INVARIANTMOMENTS_HH_

#include <vector>
#include <string>
#include <cassert>

#include "cnpy.h"

#include "MomentPreprocessor.hh"
#include "InvariantMomentCalculator.hh"

#include "ParticleDistances.hh"

template <class Event>
class InvariantMoments : public frw::AbsFrameworkAnalyzer<Event>
{
private:
	typedef typename Event::clust_seq_type::cluster_type cluster_type;

public:
	typedef Event event_type;
	typedef frw::AbsFrameworkAnalyzer<Event> Base;
	
	inline InvariantMoments(const std::string& label)
		: Base(label)
	{
	}
	
	inline virtual ~InvariantMoments() override {}
	
	inline virtual InvariantMoments* clone() const override
	{
		return new InvariantMoments(*this);
	}

private:
	
	void getLeaves(int node,
			const std::vector<cluster_type>& history,
			std::vector<int>& leaves)
	{
		if (node < 0 || node >= static_cast<int>(history.size()))
			return;
		
		const auto& c = history[node];
		
		int p1 = c.parent1();
		int p2 = c.parent2();
		
		if (p1 < 0 && p2 < 0)
		{
			leaves.push_back(node);
			return;
		}
		
		if (p1 >= 0)
			getLeaves(p1, history, leaves);
		
		if (p2 >= 0)
			getLeaves(p2, history, leaves);
	}

public:
	
	inline virtual bool analyze(const Event& evt) override
	{
		assert(evt.diffusionSequenceReady);
		assert(evt.simpleDiffusionJetsReady);
		
		const auto& history = evt.diffusionSequence.clustHist();
		const unsigned size = history.size();
		
		std::string evtnum = std::to_string(evt.number());
		
		std::vector<double> im;
		std::vector<double> imdist;
		
		constexpr int NMAX = 5;
		constexpr int MMAX = 5;
		
		for (unsigned i = 0; i < size; ++i)
		{
			// Collect jet constituents
			std::vector<int> jetParticles;
			getLeaves(i, history, jetParticles);
			
			if (jetParticles.size() < 2)
				continue;
			// Convert to generic points
			std::vector<RawMomentPoint> particles;
			particles.reserve(jetParticles.size());
			
			for (int id : jetParticles)
			{
				const auto& p = history[id].p();
				particles.push_back({p.eta(), p.phi(), p.pt()});
			}
			
			std::vector<MomentPoint> normalized = MomentPreprocessor::preprocess(particles);
			
			Moments moments = InvariantMomentCalculator::compute(normalized, NMAX, MMAX);
			
			Invariants invariants = InvariantMomentCalculator::computeInvariants(moments);
			
			for (double invariant : invariants.value)
			{
				im.push_back(invariant);
				imdist.push_back(history[i].dist());
			}
		}
			
		if (!im.empty())
		{
			cnpy::npy_save("npyarrays/im/" + evtnum + "_im.npy", &im[0], {im.size()}, "w");
		}
			
		if (!imdist.empty())
		{
			cnpy::npy_save("npyarrays/imdist/" + evtnum + "_imdist.npy", &imdist[0], {imdist.size()}, "w");
		}
			
		return true;
	}
};

#endif // INVARIENTMOMENTS_HH_
