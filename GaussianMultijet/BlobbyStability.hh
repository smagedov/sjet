#ifndef BLOBBYSTABILITY_HH_
#define BLOBBYSTABILITY_HH_

#include <cmath>
#include <string>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "cnpy.h"
#include "ParticleDistances.hh"

template <class Event>
class BlobbyStability : public frw::AbsFrameworkAnalyzer<Event>
{
public:
        typedef Event event_type;
        typedef frw::AbsFrameworkAnalyzer<Event> Base;

        inline BlobbyStability(const std::string& i_label)
                : Base(i_label)
        {
        }

        inline virtual ~BlobbyStability() override {}

        inline virtual BlobbyStability* clone() const override
                {return new BlobbyStability(*this);}

        inline virtual bool analyze(const Event& evt) override
        {
                assert(evt.diffusionSequenceReady);
                assert(evt.simpleDiffusionJetsReady);

                const unsigned nJetNodes = evt.simpleDiffusionNodes.size();
                //const unsigned nJets = evt.genJets.size();

                std::string evtnum = std::to_string(evt.number());
                std::vector<unsigned> nParts;
                unsigned curParts = 0;
                const unsigned nGenClus = evt.genClusters.size();
                for (unsigned i=0; i<nGenClus; ++i) {
                        if (!evt.genClusters[i].empty()) {
                                curParts = curParts + evt.genClusters[i].size();
                                nParts.push_back(curParts);
                        }
                }

                assert(nJetNodes == evt.simpleDiffusionJets.size());
                const std::vector<cluster_type>& history = evt.diffusionSequence.clustHist();
                const unsigned size = history.size();

		std::vector<double> stabdists;

                for (unsigned i=0; i<size; ++i) {
                        stabdists.push_back(history[i].dist());
                }

                cnpy::npy_save("npyarrays/stabdists/" + evtnum + "_stabdists.npy", &stabdists[0], {stabdists.size()}, "w");

		//Adjecency List
		std::vector<std::vector<int>> adj(size);
		std::vector<std::unordered_set<int>> tmp(size);
		for (unsigned i = 0; i < size; ++i) {
			int p1 = history[i].parent1();
			int p2 = history[i].parent2();
			if (p1 >= 0 && p1 < (int)size) {
				tmp[i].insert(p1);
				tmp[p1].insert(i);
			}
			if (p2 >= 0 && p2 < (int)size) {
				tmp[i].insert(p2);
				tmp[p2].insert(i);
			}
		}
		for (unsigned i = 0; i < size; ++i) adj[i] = std::vector<int>(tmp[i].begin(), tmp[i].end());

                //Defining Initial Field
                std::vector<double> F0(size, 0.0);
                for (unsigned i = 0; i < size; ++i) {
                        int parent1 = history[i].parent1();
                        int parent2 = history[i].parent2();

                        if (parent1 < 0 || parent2 < 0) continue;

                        int parent = (history[parent1].p().pt() >= history[parent2].p().pt()) ? parent1 : parent2;

                        double m0  = history[parent].p().m();
                        double pt0 = history[parent].p().pt();

                        if (m0 == 0 || pt0 == 0) continue;
                        double dm = std::log((history[i].p().m() + 1e-12) / (m0 + 1e-12));
			double dpt = std::log((history[i].p().pt() + 1e-12) / (pt0 + 1e-12));
			double inst = std::tanh(dm + dpt);

                        F0[i] = inst;
                }

		//Define Diffusion Step
		auto diffuse_step = [&](const std::vector<double>& in, std::vector<double>& out, double eps) {
			for (unsigned i = 0; i < size; ++i) {
				int deg = adj[i].size();
				if (deg == 0) {
					out[i] = in[i];
					continue;
				}
				double lap = 0.0;
				for (int j : adj[i]) lap += (in[j] - in[i]);
				out[i] = in[i];
				for (int j : adj[i]) out[i] += (in[j] - in[i]) / deg;
				out[i] += eps * lap;
			}
		};

		std::vector<double> scales = {0.1, 0.3, 0.6, 1.0};

		//Unified Wavelet + Heat Response
		std::vector<std::vector<double>> Fscale(scales.size(), std::vector<double>(size));
		for (size_t s = 0; s < scales.size(); ++s) {
			double eps = 0.05;
			int nSteps = (int)(scales[s] / eps);
			std::vector<double> F = F0;
			std::vector<double> next(size);
			
			for (int t = 0; t < nSteps; ++t) {
				diffuse_step(F, next, eps);
				F = next;
			}
			Fscale[s] = F;
		}

		//Unified Wavelet Operator
		std::vector<std::vector<double>> W(scales.size()-1, std::vector<double>(size));
		for (size_t s = 0; s + 1 < scales.size(); ++s) {
			for (unsigned i = 0; i < size; ++i) {
				W[s][i] = Fscale[s][i] - Fscale[s+1][i];
			}
		}

		//Blobbiness
		std::vector<double> stability(size, 0.0);
		for (unsigned i = 0; i < size; ++i) {
			double max_wavelet = 0.0;
			for (size_t s = 0; s < W.size(); ++s) max_wavelet = std::max(max_wavelet, std::abs(W[s][i]));
			double support = std::abs(Fscale.back()[i]);
			stability[i] = max_wavelet * support;
		}

                cnpy::npy_save("npyarrays/stability/" + evtnum + "_stability.npy", &stability[0], {stability.size()}, "w");

		stability.clear();

		return true;
	}

private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;

};

#endif // BLOBBYSTABILITY_HH_

