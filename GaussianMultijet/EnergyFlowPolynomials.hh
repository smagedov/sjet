#ifndef ENERGYFLOWPOLYNOMIALS_HH_
#define ENERGYFLOWPOLYNOMIALS_HH_

#include <cmath>
#include <string>

#include "cnpy.h"
#include "ParticleDistances.hh"

struct Edge {
    int a;
    int b;
};

struct Graph {
    int V;
    std::vector<Edge> edges;
};

template <class Event>
class EnergyFlowPolynomials : public frw::AbsFrameworkAnalyzer<Event>
{
private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;
public:
        typedef Event event_type;
        typedef frw::AbsFrameworkAnalyzer<Event> Base;

        inline EnergyFlowPolynomials(const std::string& i_label)
                : Base(i_label)
        {
        }

        inline virtual ~EnergyFlowPolynomials() override {}

        inline virtual EnergyFlowPolynomials* clone() const override
                {return new EnergyFlowPolynomials(*this);}

	void getLeaves(int node, 
			const std::vector<cluster_type>& history,
			std::vector<int>& leaves)
	{
		if (node < 0 || node >= static_cast<int>(history.size())) {
			return;
		}
		
		const auto& c = history[node];
		
		int p1 = c.parent1();
		int p2 = c.parent2();
		
		if (p1 < 0 && p2 < 0) {
			leaves.push_back(node);
			return;
		}
		
		if (p1 >= 0) {
			getLeaves(p1, history, leaves);
		}
		
		if (p2 >= 0) {
			getLeaves(p2, history, leaves);
		}
	}

	std::vector<double> computeZ(const std::vector<int>& jet_particles, const std::vector<cluster_type>& history)
	{
		double Etot = 0.0;
		
		for (int id : jet_particles) {
			Etot += history[id].p().e();
		}
		
		std::vector<double> z;
		z.reserve(jet_particles.size());
		
		for (int id : jet_particles) {
			z.push_back(history[id].p().e()/Etot);
		}
		
		return z;
	}

	std::vector<std::vector<double>> computeDRMatrix(const std::vector<int>& jet_particles, const std::vector<cluster_type>& history)
	{
		size_t N = jet_particles.size();
		std::vector<std::vector<double>> dr(N, std::vector<double>(N));
		const ParticleDeltaR dRcalculator;
		
		for (size_t i = 0; i < N; i++) {
			for (size_t j = 0; j < N; j++) {
				dr[i][j] = dRcalculator(history[jet_particles[i]].p(),history[jet_particles[j]].p());
			}
		}
		
		return dr;
	}

	double evaluateEFP(const std::vector<int>& jet_particles, const Graph& G, const std::vector<cluster_type>& history)
	{
		auto z  = computeZ(jet_particles, history);
		auto dr = computeDRMatrix(jet_particles, history);
		
		int N = jet_particles.size();
		
		std::vector<int> idx(G.V);
		double result = 0.0;
		
		std::function<void(int)> dfs = [&](int v)
		{
			if (v == G.V)
			
			{
				double term = 1.0;
				// energy terms
				for (int a = 0; a < G.V; a++) {
					term *= z[idx[a]];
				}
				
				// angular terms
				for (auto& e : G.edges) {
					term *= dr[idx[e.a]][idx[e.b]];
				}
				
				result += term;
				return;
			}
			
			for (int i = 0; i < N; i++) {
				idx[v] = i;
				dfs(v + 1);
			}
		};
		
		dfs(0);
		return result;
	}


        inline virtual bool analyze(const Event& evt) override
        {
                assert(evt.diffusionSequenceReady);
                assert(evt.simpleDiffusionJetsReady);

                const unsigned nJetNodes = evt.simpleDiffusionNodes.size();
                //const unsigned nJets = evt.genJets.size();

                std::string evtnum = std::to_string(evt.number());
                Graph G;
		G.V = 2;
		G.edges = {{0,1}};

                assert(nJetNodes == evt.simpleDiffusionJets.size());
                const std::vector<cluster_type>& history = evt.diffusionSequence.clustHist();
                const unsigned size = history.size();

		std::vector<double> efp;
		std::vector<double> efpdist;
		for (unsigned i=0; i<size; ++i) {
			std::vector<int> jet_particles;
			getLeaves(i, history, jet_particles);
			if (jet_particles.size() < 2) {
				continue;
			}

			double val = evaluateEFP(jet_particles, G, history);
			efp.push_back(val);
			efpdist.push_back(history[i].dist());
		}

                cnpy::npy_save("npyarrays/efp/" + evtnum + "_efp.npy", &efp[0], {efp.size()}, "w");
                cnpy::npy_save("npyarrays/efpdist/" + evtnum + "_efpdist.npy", &efpdist[0], {efpdist.size()}, "w");

		efp.clear();
		efpdist.clear();

		return true;
	}
};

#endif // ENERGYFLOWPOLYNOMIALS_HH_


