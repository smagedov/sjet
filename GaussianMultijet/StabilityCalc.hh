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
		std::vector<double> stability(size);
		for (unsigned i=0; i<nParts[nGenClus-1]; ++i) {
			stability[i] = 0.0;
		}

                for (unsigned i=0; i<size; ++i) {
                        stabdists.push_back(history[i].dist());
                }

                cnpy::npy_save("npyarrays/" + evtnum + "_stabdists.npy", &stabdists[0], {stabdists.size()}, "w");

		//Stability going up:
		//const ParticleDeltaR dRcalculator;
		for (unsigned i=nParts[nGenClus-1]; i<size; ++i) {
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
				double distratio = history[parent].dist()/history[i].dist();
				//std::cout << i << " " << parent << " " << distratio << std::endl;
				//double massratio = history[parent].p().m()/history[i].p().m();
				double ptratio = history[parent].p().pt()/history[i].p().pt();
				//double deltar = dRcalculator(history[parent].p(), history[i].p());
				stab = stability[parent] + ptratio*log(distratio);
			}
			stability[i] = stab;
		}

		//Stability going down:
		for (unsigned i=size-1; i>=nParts[nGenClus-1]; --i) {
			int daughter = history[i].daughter();
			double stab = 0.0;
			if (history[i].dist() != -1) {
				double distratio = history[i].dist()/history[daughter].dist();
				double ptratio = history[i].p().pt()/history[daughter].p().pt();
				stab = stability[daughter] + ptratio*log(distratio);
			}
			stability[i] += stab;
		}

                cnpy::npy_save("npyarrays/" + evtnum + "_stability.npy", &stability[0], {stability.size()}, "w");


		return true;
	}

private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;

};

#endif // STABILITYCALC_HH_

