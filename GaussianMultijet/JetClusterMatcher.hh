#ifndef JETCLUSTERMATCHER_HH_
#define JETCLUSTERMATCHER_HH_

#include <cmath>

#include "cnpy.h"
#include "ParticleDistances.hh"

template <class Event>
class JetCluserMatcher : public frw::AbsFrameworkAnalyzer<Event>
{
public:
        typedef Event event_type;
        typedef frw::AbsFrameworkAnalyzer<Event> Base;

        inline JetClusterMatcher(const std::string& i_label)
                : Base(i_label)
        {
        }

        inline virtual ~JetClusterMatcher() override {}

        inline virtual JetClusterMatcher* clone() const override
                {return new JetClusterMatcher(*this);}

        inline virtual bool analyze(const Event& evt) override
        {
                assert(evt.diffusionSequenceReady);
                assert(evt.simpleDiffusionJetsReady);

                const unsigned nJetNodes = evt.simpleDiffusionNodes.size();
                const unsigned nJets = evt.genJets.size();

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
		std::vector<int> mask(size);

		std::vector<unsigned> closestjeti;
                for (unsigned i=0; i<nJets; ++i) {
                        unsigned closesti = 0;
                        double closestdist = dRcalculator(evt.genJets[i], history[nParts[nJets-1]].p());
                        for (unsigned j=nParts[nJets]; j<size; ++j) {
				if (mask[j] != 0) {
                                	double tmp = dRcalculator(evt.genJets[i], history[j].p());
                                	if (tmp < closestdist) {
                                        	closestdist = tmp;
                                        	closesti = j;
                                	}
				}
                        }
                        closestjeti.push_back(closesti);
			std::cout << "parent1 id:" << history[closesti].parent1() << "parent2 id:" << history[closesti].parent2() << std::endl;

                }

                return true;
        }

private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;

};

#endif // JETHISTORYCOPY_HH_

