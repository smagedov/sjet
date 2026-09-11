#ifndef INVARIANTMOMENTS_HH_
#define INVARIANTMOMENTS_HH_

#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include <utility>

#include "cnpy.h"

#include "MomentAccumulator.hh"
#include "getInvariants.hh"

#include "ParticleDistances.hh"
#include "GlobalJetMatcher.hh"

namespace {
	struct RawMomentPoint
	{
		double eta;
		double phi;
		double weight;
	};


        class FCN
        {
                public:
                        inline double operator()(const double /*x*/) const
                        {
                                return 1.0;
                        }
        };

	bool approximatelyEqual(double a, double b)
	{
		const double abs_eps = 1e-8;
		const double rel_eps = 1e-6;
		
		if (std::abs(a-b) < abs_eps)
			return true;

		double denom = std::max(std::abs(a), std::abs(b));
		
		return std::abs(a-b)/denom < rel_eps;
	}

}

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
	void disableAncestors(int idx, 
			const std::vector<cluster_type>& history,
			std::vector<bool>& available)
	{
		if (idx < 0 || !available[idx])
			return;
		
		available[idx] = false;
		
		disableAncestors(history[idx].parent1(), history, available);
		disableAncestors(history[idx].parent2(), history, available);
	}
	
	void disableDaughters(int idx,
			const std::vector<cluster_type>& history,
			std::vector<bool>& available)
	{
		int daughter = history[idx].daughter();
		
		if (daughter >= static_cast<int>(history.size()) || !available[daughter])
			return;
		
		available[daughter] = false;
		disableDaughters(daughter, history, available);
	}

	void collectHistoryIndices(int idx,
                           const std::vector<cluster_type>& history,
                           std::vector<int>& indices,
                           std::vector<bool>& visited)
	{
		if (idx < 0 || idx >= static_cast<int>(history.size()))
			return;
		
		if (visited[idx])
			return;
		
		visited[idx] = true;
		indices.push_back(idx);
		
		int p1 = history[idx].parent1();
		int p2 = history[idx].parent2();
		
		int bestParent = -1;
		
		if (p1 >= 0 && p2 >= 0)
		{
			bestParent = (history[p1].p().pt() > history[p2].p().pt()) ? p1 : p2;
		}
		else if (p1 >= 0)
		{
			bestParent = p1;
		}
		else if (p2 >= 0)
		{
			bestParent = p2;
		}
		
		// Follow the leading-pT parent branch
		collectHistoryIndices(bestParent, history, indices, visited);
		
		// Then follow the daughter
		collectHistoryIndices(history[idx].daughter(), history, indices, visited);
	}

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

	std::vector<double> calculateInvariants(const std::vector<RawMomentPoint>& particles, double sigma)
	{
		long double eta_acc = 0.0L;
		long double sin_acc = 0.0L;
		long double cos_acc = 0.0L;
		long double w_acc = 0.0L;
		
		for (const auto& p : particles) {
			eta_acc += p.eta * p.weight;
			sin_acc += std::sin(p.phi) * p.weight;
			cos_acc += std::cos(p.phi) * p.weight;
			w_acc += p.weight;
		}
		
		if (w_acc == 0.0L)
			return {};
		
		double eta_center = eta_acc / w_acc;
		long double phi_center_sin = sin_acc;
		long double phi_center_cos = cos_acc;
		
		double phi_center = atan2(static_cast<double>(phi_center_sin), static_cast<double>(phi_center_cos));
		
		MomentAccumulator<FCN> acc(FCN(),
				eta_center,
				0.0,
				sigma);
		
		for (const auto& p : particles) {
			double dphi = p.phi - phi_center;
			while(dphi > M_PI)
				dphi -= 2.0*M_PI;
			while(dphi <= -M_PI)
				dphi += 2.0*M_PI;

			acc.accumulate(p.eta,
					dphi,
					p.weight);
		}
		
		constexpr unsigned DIM = ghm::MAXDEGP1;
		double moments[DIM][DIM];
		acc.getMoments(moments);
		
		std::vector<double> invs;
		for (unsigned i = 1; i <= 18; ++i)
			invs.push_back(getInvariants(i,moments));
		return invs;
	}

public:
	inline virtual bool analyze(const Event& evt) override
	{
		assert(evt.diffusionSequenceReady);
		assert(evt.simpleDiffusionJetsReady);
		
		const auto& history = evt.diffusionSequence.clustHist();
		const unsigned size = history.size();
		
		const unsigned nJets = evt.genJets.size();
		
		std::string evtnum = std::to_string(evt.number());
		
		// ------------------------------------------------------------
		// Match each oracle/gen jet to its closest available history node.
		// This is the same matching procedure used in JetClusterMatcher.
		// ------------------------------------------------------------

		const ParticleDeltaR dRcalculator;
		
		std::vector<bool> available(size, true);
		std::vector<unsigned> closestjeti;
		
		const double alpha = 1.0;
		const double beta = 1.0;
		const unsigned maxCandidates = size;
		
		auto result = GlobalJetMatching::match(evt.genJets,
				history,
				nJets,
				alpha,
				beta,
				dRcalculator,
				// This tells the matcher which clusters are incompatible
				// with a selected cluster.
				[&](unsigned j, std::vector<bool>& available)
				{
				disableAncestors(j, history, available);
				disableDaughters(j, history, available);
				},
				maxCandidates,
				// Print top candidate clusters for debugging.
				true
				);

		closestjeti.clear();
		closestjeti.resize(nJets, size);
		
		if (!result.found) {
			std::cout << "ERROR: Could not find globally compatible " << nJets << "-jet matching using " << maxCandidates << " candidates per oracle jet." << std::endl;
		} else {
			std::cout << "\nGLOBAL MATCH TOTAL SCORE: " << result.totalScore << std::endl;
			bool goodGlobalMatch = true;
			const double maxDeltaR = 0.30;
			const double maxPtRatio = 0.30;
			const double maxMassRatio = 0.30;
			for (unsigned i = 0; i < nJets; ++i) {
				const unsigned j = result.clusterForOracle[i];
				const double deltaR = dRcalculator(evt.genJets[i], history[j].p());
				const double ptRatio = std::abs(std::log(evt.genJets[i].pt() / history[j].p().pt()));
				const double massRatio = std::abs(std::log((evt.genJets[i].m() + 1) / (history[j].p().m() + 1)));
				const double score = deltaR + alpha*ptRatio + beta*massRatio;
				
				const bool passed = deltaR < maxDeltaR && ptRatio < maxPtRatio && massRatio < maxMassRatio;
				if (!passed) {
					goodGlobalMatch = false;
				}
				std::cout 
					<< "Oracle Jet: " << i 
					<< " closest id: " << j  
					<< " max size: " << size 
					<< " score: " << score  
					<< " deltaR: " << deltaR  
					<< " PT Ratio: " << ptRatio 
					<< " mass Ratio: " << massRatio 
					<< std::endl;
			}

                        //if (!goodGlobalMatch) {
			//	std::cout << "GLOBAL MATCH REJECTED: " << "at least one oracle jet failed quality cuts." << std::endl;
			//	return true;
			//}

			for (unsigned i = 0; i < nJets; ++i) {
				closestjeti[i] = result.clusterForOracle[i];
			}

		}
		
		// ------------------------------------------------------------
		// Process each matched jet independently.
		// ------------------------------------------------------------

		std::vector<double> debug_deltaR;
                std::vector<double> debug_ptRatio;
                std::vector<double> debug_score;
                std::vector<double> debug_oraclePt;
                std::vector<double> debug_clusterPt;
		std::vector<double> debug_dist;
		for (unsigned jetIndex = 0; jetIndex < closestjeti.size(); ++jetIndex)
		{
			const int jetNode = closestjeti[jetIndex];
			double mdr = dRcalculator(evt.genJets[jetIndex], history[jetNode].p());
			double mptrat = evt.genJets[jetIndex].pt()/history[jetNode].p().pt();
			double score = mdr + std::abs(std::log(mptrat));
			double orpt = evt.genJets[jetIndex].pt();
			double mpt = history[jetNode].p().pt();
			double dist = history[jetNode].dist();

			debug_deltaR.push_back(mdr);
			debug_ptRatio.push_back(mptrat);
			debug_score.push_back(score);
			debug_oraclePt.push_back(orpt);
			debug_clusterPt.push_back(mpt);
			debug_dist.push_back(dist);


			//if ((mdr < 0.2) && ((mptrat < 1.1) && (mptrat > 0.9))) {
			//Obtain a vector of oracle and matched imaginary moments
			std::vector<double> oracleIm;
			std::vector<double> matchedIm;
				
			std::vector<int> oracleParticles;
			std::vector<int> matchedParticles;
			getLeaves(jetNode, history, matchedParticles);
			std::vector<RawMomentPoint> orparticles;
			std::vector<RawMomentPoint> maparticles;
			for (int id : evt.genClusters[jetIndex])
			{
				const Pythia8::Particle& p = (*evt.pythiaEvent)[id];
				rk::P4 part = rk::P4(p.pT()*geom3::Vector3(cos(p.phi()), sin(p.phi()), sinh(p.eta())), p.m());
			
				orparticles.push_back({part.eta(), part.phi(), part.pt()});
			}
	
			for (int id : matchedParticles)
			{
       	                       	const auto& p = history[id].p();
		
				maparticles.push_back({p.eta(), p.phi(), p.pt()});
	                }
	                double alpha = 1.0;
	                double sigma = history[jetNode].dist();
	
       	                // Calculate all 18 invariants.
       	                std::vector<double> oracleinvs = calculateInvariants(orparticles, alpha*sigma);
			std::vector<double> matchedinvs = calculateInvariants(maparticles, alpha*sigma);
	
			// Store the 18 invariants consecutively.
       	                for (unsigned i = 0; i < oracleinvs.size(); ++i)
			{
       	                        oracleIm.push_back(oracleinvs[i]);
				matchedIm.push_back(matchedinvs[i]);
			}
	
			std::cout << std::to_string(jetIndex) << " Jet Passed: We good." << std::endl;
       	                if (!oracleIm.empty())
       	                {
       	                       	cnpy::npy_save("npyarrays/debugging/im_oracle/" + evtnum + "_jet" + std::to_string(jetIndex) + "_im.npy", &oracleIm[0], {oracleIm.size()}, "w");
       	                }
       	                if (!oracleIm.empty())
       	                {
       	                        cnpy::npy_save("npyarrays/debugging/im_matched/" + evtnum + "_jet" + std::to_string(jetIndex) + "_im.npy", &matchedIm[0], {matchedIm.size()}, "w");
			}

			
			// Get the ordered history of this particular jet.
			std::vector<int> historyIndices;
			std::vector<bool> visited(size, false);
			
			collectHistoryIndices(jetNode, history, historyIndices, visited);
			
			// One row per history step.
			// Each row contains the 18 invariant moments.
			std::vector<double> jetIm;
			std::vector<double> jetImDist;
			std::vector<double> jetRadii;
			
			jetIm.reserve(historyIndices.size() * 18);
			jetImDist.reserve(historyIndices.size());
			jetRadii.reserve(historyIndices.size());
			
			// --------------------------------------------------------
			// Calculate invariants for every node in this jet's history.
			// --------------------------------------------------------

			for (int node : historyIndices)
			{
				// Get all leaf particles belonging to this history node.
				std::vector<int> jetParticles;
				getLeaves(node, history, jetParticles);
				
				// Need at least two particles.
				if (jetParticles.size() < 2)
					continue;
				
				// Convert leaves to RawMomentPoints.
				std::vector<RawMomentPoint> particles;
				particles.reserve(jetParticles.size());
				
				for (int id : jetParticles)
				{
					const auto& p = history[id].p();
					
					particles.push_back({p.eta(), p.phi(), p.pt()});
				}
		
				// ------------------------------------------------------------
                		// Calculate pT-weighted 68% containment radius.
                		// ------------------------------------------------------------
				
			        long double eta_acc = 0.0L;
				long double sin_acc = 0.0L;
				long double cos_acc = 0.0L;
                		long double w_acc = 0.0L;

                		for (const auto& p : particles) {
                        		eta_acc += p.eta * p.weight;
                        		sin_acc += std::sin(p.phi) * p.weight;
                        		cos_acc += std::cos(p.phi) * p.weight;
                        		w_acc += p.weight;
                		}

                		if (w_acc == 0.0L)
                        		return {};

                		double eta_center = eta_acc / w_acc;
                		long double phi_center_sin = sin_acc;
                		long double phi_center_cos = cos_acc;

                		double phi_center = atan2(static_cast<double>(phi_center_sin), static_cast<double>(phi_center_cos));

                		std::vector<std::pair<double, double>> radialPoints;
                		radialPoints.reserve(particles.size());

                		for (const auto& p : particles) {
                        		double deta = p.eta - eta_center;
                        		double dphi = p.phi - phi_center;

                        		while (dphi > M_PI)
                                		dphi -= 2.0 * M_PI;
                        		while (dphi <= -M_PI)
                                		dphi += 2.0 * M_PI;

                        		double radius = std::sqrt(deta * deta + dphi * dphi);
                        		radialPoints.emplace_back(radius, p.weight);
                		}
                		// Sort from smallest radius to largest radius.
                		std::sort(radialPoints.begin(), radialPoints.end(), [](const auto& a, const auto& b) {
                                		return a.first < b.first;
                                		});

                		// Find the radius containing 86.5% of the total pT.
                		const long double targetWeight = 0.865L * w_acc;

                		long double accumulatedWeight = 0.0L;
                		double radii = 0.0;

                		for (const auto& point : radialPoints) {
                        		accumulatedWeight += point.second;
                        		radii = point.first;
                        		if (accumulatedWeight >= targetWeight)
                                		break;
                		}


				double alpha = 1.0;
				double sigma = history[node].dist();
				
				// Calculate all 18 invariants.
				std::vector<double> invs = calculateInvariants(particles, alpha*sigma);
				
				// Store the 18 invariants consecutively.
				for (double inv : invs)
					jetIm.push_back(inv);
				
				// Distance corresponding to this history step.
				jetImDist.push_back(history[node].dist());

				// Characteristic distance.
				jetRadii.push_back(radii);
			}
			
			// --------------------------------------------------------
			// Save this jet's history separately.
			//
			// jetIm shape:  [number_of_history_steps * 18]
			//
			// In Python it can be reshaped to:
			//     (number_of_history_steps, 18)
			// --------------------------------------------------------

			if (!jetIm.empty())
			{
				cnpy::npy_save("npyarrays/im/" + evtnum + "_jet" + std::to_string(jetIndex) + "_im.npy", &jetIm[0], {jetIm.size()}, "w");
			}
			
			if (!jetImDist.empty())
			{
				cnpy::npy_save("npyarrays/distcomp/imdist/" + evtnum + "_jet" + std::to_string(jetIndex) + "_imdist.npy", &jetImDist[0], {jetImDist.size()}, "w");
			}

                        if (!jetImDist.empty())
                        {
                                cnpy::npy_save("npyarrays/distcomp/radii/" + evtnum + "_jet" + std::to_string(jetIndex) + "_radii.npy", &jetRadii[0], {jetRadii.size()}, "w");
                        }

		}

		std::vector<double> invmass;
		std::vector<double> recmass;
		std::vector<double> clusmass;
		for (int i=0; i<evt.invMasses.size(); ++i) {
			invmass.push_back(evt.invMasses[i]);
			recmass.push_back(evt.recMasses[i]);
		}

		const rk::P4& b     = history[closestjeti[0]].p();
		const rk::P4& wq1   = history[closestjeti[1]].p();
		const rk::P4& wq2   = history[closestjeti[2]].p();
		const rk::P4& bbar  = history[closestjeti[3]].p();
		const rk::P4& wq3   = history[closestjeti[4]].p();
		const rk::P4& wq4   = history[closestjeti[5]].p();
			
		const rk::P4 wPlus  = wq1 + wq2;
		const rk::P4 wMinus = wq3 + wq4;		
		const rk::P4 top     = b    + wPlus;
		const rk::P4 antiTop = bbar + wMinus;

		clusmass.push_back(wPlus.m());
		clusmass.push_back(wMinus.m());
		clusmass.push_back(top.m());
		clusmass.push_back(antiTop.m());

		cnpy::npy_save("npyarrays/debugging/invmass/" + evtnum + "_invmass.npy", &invmass[0], {invmass.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/recmass/" + evtnum + "_recmass.npy", &recmass[0], {recmass.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/clusmass/" + evtnum + "_clusmass.npy", &clusmass[0], {clusmass.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/deltar/" + evtnum + "_deltar.npy", &debug_deltaR[0], {debug_deltaR.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/ptratio/" + evtnum + "_ptratio.npy", &debug_ptRatio[0], {debug_ptRatio.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/oraclept/" + evtnum + "_oraclept.npy", &debug_oraclePt[0], {debug_oraclePt.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/matchedpt/" + evtnum + "_matchedpt.npy", &debug_clusterPt[0], {debug_clusterPt.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/dist/" + evtnum + "_dist.npy", &debug_dist[0], {debug_dist.size()}, "w");
		cnpy::npy_save("npyarrays/debugging/score/" + evtnum + "_score.npy", &debug_score[0], {debug_score.size()}, "w");
		
		return true;
	}
};

#endif // INVARIENTMOMENTS_HH_
