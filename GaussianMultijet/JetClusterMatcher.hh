#ifndef JETCLUSTERMATCHER_HH_
#define JETCLUSTERMATCHER_HH_

#include <cmath>

#include "cnpy.h"
#include "ParticleDistances.hh"

template <class Event>
class JetClusterMatcher : public frw::AbsFrameworkAnalyzer<Event>
{
private:
        typedef typename Event::clust_seq_type::cluster_type cluster_type;
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

	void disableAncestors(int idx, const std::vector<cluster_type>& history, std::vector<bool>& available)
	{
		if (idx < 0 || !available[idx]) {
			return;
		}

		available[idx] = false;
		disableAncestors(history[idx].parent1(), history, available);
		disableAncestors(history[idx].parent2(), history, available);
	}

        void disableDaughters(int idx, const std::vector<cluster_type>& history, std::vector<bool>& available)
        {
		int daughter = history[idx].daughter();
                if (daughter >= history.size() || !available[daughter]) {
                        return;
		}
		
                available[daughter] = false;
                disableAncestors(daughter, history, available);
        }

	void collectHistoryPts(int idx, const std::vector<cluster_type>& history, 
			std::vector<float>& pts, std::vector<float>& ms, 
			std::vector<float>& dists, std::vector<bool>& visited)
	{
		    if (idx < 0 || idx >= history.size()) return;
		    if (visited[idx]) return;

		    visited[idx] = true;

		    pts.push_back(history[idx].p().pt());
		    ms.push_back(history[idx].p().m());
		    dists.push_back(history[idx].dist());

		    collectHistoryPts(history[idx].parent1(), history, pts, ms, dists, visited);
		    collectHistoryPts(history[idx].parent2(), history, pts, ms, dists, visited);
		    collectHistoryPts(history[idx].daughter(), history, pts, ms, dists, visited);;
	}

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

		std::cout << "Total Oracle Jets: " << nJets << " History Size: " << size << std::endl;
		const ParticleDeltaR dRcalculator;
		std::vector<bool> available(size, true);
		std::vector<unsigned> closestjeti;
		double alpha = 0.04;
		for (unsigned i=0; i<nJets; ++i) {
                        unsigned closesti = 0;
                        double closestdist = dRcalculator(evt.genJets[i], history[0].p()) + alpha*std::log(evt.genJets[i].pt()/history[0].p().pt());
                        for (unsigned j=1; j<size; ++j) {
				if (available[j]) {
                                	double tmp = dRcalculator(evt.genJets[i], history[j].p()) + alpha*std::log(evt.genJets[i].pt()/history[j].p().pt());
                                	if (tmp < closestdist) {
                                        	closestdist = tmp;
                                        	closesti = j;
                                	}
				}
                        }
			disableAncestors(closesti, history, available);
			disableDaughters(closesti, history, available);

                        closestjeti.push_back(closesti);
			std::cout << "Oracle Jet: " << i << " closest id:" << closesti << " closest dist:" << closestdist << std::endl;

                }

		std::vector<std::vector<float>> jetHistoryPts;
		std::vector<std::vector<float>> jetHistoryMs;
		std::vector<std::vector<float>> jetHistoryDists;
		for (unsigned jeti : closestjeti) {
			std::vector<float> pts;
			std::vector<float> ms;
			std::vector<float> dists;
			std::vector<bool> visited(history.size(), false);

			collectHistoryPts(jeti, history, pts, ms, dists, visited);
			jetHistoryPts.push_back(pts);
			jetHistoryMs.push_back(ms);
			jetHistoryDists.push_back(dists);
		}

		std::vector<float> flatPts;
		std::vector<float> flatMs;
		std::vector<float> flatDists;
		std::vector<int> offsets;
		
		int cursor = 0;
		
		for (size_t i=0; i<jetHistoryPts.size(); ++i) {
			const auto& ptsRow = jetHistoryPts[i];
			const auto& msRow = jetHistoryMs[i];
			const auto& distRow = jetHistoryDists[i];
			offsets.push_back(cursor);
			flatPts.insert(flatPts.end(), ptsRow.begin(), ptsRow.end());
			flatMs.insert(flatMs.end(), msRow.begin(), msRow.end());
			flatDists.insert(flatDists.end(), distRow.begin(), distRow.end());
			cursor += ptsRow.size();
		}
		offsets.push_back(cursor);

		std::string evtnum = std::to_string(evt.number());
		cnpy::npy_save("npyarrays/jethistory/flatDists/" + evtnum + "_flatDists.npy", &flatDists[0], {flatDists.size()}, "w");
		cnpy::npy_save("npyarrays/jethistory/flatPts/" + evtnum + "_flatPts.npy", &flatPts[0], {flatPts.size()}, "w");
		cnpy::npy_save("npyarrays/jethistory/flatMs/" + evtnum + "_flatMs.npy", &flatMs[0], {flatMs.size()}, "w");
		cnpy::npy_save("npyarrays/jethistory/offsets/" + evtnum + "_offsets.npy", &offsets[0], {offsets.size()}, "w");

		jetHistoryPts.clear();
		flatPts.clear();
		offsets.clear();
                closestjeti.clear();
		return true;
        }
};
#endif // JETHISTORYCOPY_HH_

