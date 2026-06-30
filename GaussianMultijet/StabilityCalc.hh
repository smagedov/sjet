#ifndef STABILITYCALC_HH_
#define STABILITYCALC_HH_

#include <cmath>
#include <string>

#include "cnpy.h"
#include "ParticleDistances.hh"

template <class Event>
class StabilityCalc : public frw::AbsFrameworkAnalyzer<Event>
{
public:
        typedef Event event_type;
        typedef frw::AbsFrameworkAnalyzer<Event> Base;

        inline StabilityCalc(const std::string& i_label)
                : Base(i_label)
        {
        }

        inline virtual ~StabilityCalc() override {}

        inline virtual StabilityCalc* clone() const override
                {return new StabilityCalc(*this);}

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
		std::vector<double> stabup(size);
		std::vector<double> stabdown(size);

                for (unsigned i=0; i<size; ++i) {
                        stabdists.push_back(history[i].dist());
                }

                cnpy::npy_save("npyarrays/stabdists/" + evtnum + "_stabdists.npy", &stabdists[0], {stabdists.size()}, "w");

		//Stability going up:
		const ParticleDeltaR dRcalculator;
		double sigma = 0.4;
		for (unsigned i=nParts[nGenClus-1]; i<size-1; ++i) {
			int parent1 = history[i].parent1();
			int parent2 = history[i].parent2();
			int parent = 0;
			double stab = 0.0;
			if (history[parent1].p().pt() >= history[parent2].p().pt()) {
				parent = parent1;
			} else {
				parent = parent2;
			}
			if (history[parent].dist() != -1) {
				double parentdist = history[parent].dist();
				double parentmass = history[parent].p().m();
				double parentpt = history[parent].p().pt();
				double partdist = history[i].dist();
				double partmass = history[i].p().m();
				double partpt = history[i].p().pt();
				double inst = (partmass - parentmass)/parentmass + (partpt - parentpt)/parentpt;
				double deltar = dRcalculator(history[parent].p(), history[i].p());
				double weight = 1 - std::exp(-(deltar*deltar)/(2*sigma*sigma));
				stab = stabup[parent] + inst*weight;
			}
			stabup[i] = stab;
		}

		//Stability going down:
		for (unsigned i=size-2; i>=nParts[nGenClus-1]; --i) {
			int daughter = history[i].daughter();
			double stab = 0.0;
			if (history[i].dist() != -1) {
                                double daughtdist = history[daughter].dist();
                                double daughtmass = history[daughter].p().m();
                                double daughtpt = history[daughter].p().pt();
                                double partdist = history[i].dist();
                                double partmass = history[i].p().m();
                                double partpt = history[i].p().pt();
                                double inst = (daughtmass - partmass)/partmass + (daughtpt - partpt)/partpt;
				double deltar = dRcalculator(history[daughter].p(), history[i].p());
				double weight = 1 - std::exp(-(deltar*deltar)/(2*sigma*sigma));
				stab = stabdown[daughter] + inst*weight;
			}
			stabdown[i] = stab;
		}

		std::vector<double>stability (size);
		for (unsigned i=0; i<size; ++i) {
			stability[i] = stabup[i] + stabdown[i];
			//std::cout << "i: " << i << " stabup: " << stabup[i] << " stabdown: " << stabdown[i] << std::endl;
		}

                cnpy::npy_save("npyarrays/stability/" + evtnum + "_stability.npy", &stability[0], {stability.size()}, "w");

		stability.clear();

		return true;
	}

private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;

};

#endif // STABILITYCALC_HH_

