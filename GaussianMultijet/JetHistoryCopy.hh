#ifndef JETHISTORYCOPY_HH_
#define JETHISTORYCOPY_HH_

#include "cnpy.h"

template <class Event>
class JetHistoryCopy : public frw::AbsFrameworkAnalyzer<Event>
{
public:
	typedef Event event_type;
	typedef frw::AbsFrameworkAnalyzer<Event> Base;

	inline JetHistoryCopy(const std::string& i_label)
		: Base(i_label)
	{
	}

	inline virtual ~JetHistoryCopy() override {}
	
	inline virtual JetHistoryCopy* clone() const override
		{return new JetHistoryCopy(*this);}
		
	inline virtual bool analyze(const Event& evt) override
	{
		assert(evt.diffusionSequenceReady);
		assert(evt.simpleDiffusionJetsReady);

		const unsigned nJetNodes = evt.simpleDiffusionNodes.size();
		const unsigned nJets = evt.genJets.size();

		std::vector<unsigned> nParts;
		unsigned curParts = 0;
		for (unsigned i=0; i<evt.genClusters.size(); ++i) {
			if (!evt.genClusters[i].empty()) {
				curParts = curParts + evt.genClusters[i].size();
				nParts.push_back(curParts);
			}
		}
		assert(nJetNodes == evt.simpleDiffusionJets.size());
		const std::vector<cluster_type>& history = evt.diffusionSequence.clustHist();

		for (unsigned i=0; i<nJets; ++i) {
			std::cout << "Jet #" << i+1 << " Total Particles: " << nParts[i] << std::endl;
		}

		std::vector<std::vector<rk::P4>> histcopy(history.size(), std::vector<rk::P4>(nJets));

		for (unsigned i=0; i<history.size(); ++i) {
			for (unsigned j=0; j<nJets; ++j) {
				if (i < nParts[j]) {
					histcopy[i][j] = history[i].p();
					break;
				}
			}
		}

		for (unsigned i=nParts[nJets-1]; i<history.size(); ++i) {
			int pind1 = history[i].parent1();
			int pind2 = history[i].parent2();
			for (unsigned j=0; j<nJets; ++j) {
				histcopy[i][j] = histcopy[pind1][j] + histcopy[pind2][j];
			}
		}

		std::vector<double> leadingPt;
		for (unsigned i = 0; i<history.size(); ++i) {
			double leadPt = 0.0;
			for (unsigned j = 0; j<nJets; ++j) {
				double tmp = histcopy[i][j].pt();
				if (tmp > leadPt) {
					leadPt = tmp;
				}
			}
			leadingPt.push_back(leadPt);
		}

		std::vector<double> dists;
		std::vector<double> ratios;
		std::vector<double> vecsumratios;
		for (unsigned i = 0; i<history.size(); ++i) {
			double vecsum = std::sqrt(history[i].p().px()*history[i].p().px() + history[i].p().py()*history[i].p().py() + history[i].p().pz()*history[i].p().pz());
			vecsumratios.push_back(leadingPt[i]/vecsum);
			ratios.push_back(leadingPt[i]/history[i].scalarPtSum());
			dists.push_back(history[i].dist());
		}
		cnpy::npy_save("npyarrays/ratios.npy", &ratios[0], {ratios.size()}, "w");
		cnpy::npy_save("npyarrays/vecsumratios.npy", &vecsumratios[0], {vecsumratios.size()}, "w");
		cnpy::npy_save("npyarrays/dists.npy", &dists[0], {dists.size()}, "w");

		for (unsigned i = 0; i<history.size(); ++i) {
			//auto& vec = evt.copySequence.clustHist()[i].p();
			//vec.insert(vec.end(), histcopy[i].begin(), histcopy[i].end());
			std::cout << "evtHist: " <<  evt.diffusionSequence.clustHist()[i].p() << " copyHist: ";
			//std::cout << "i: " << i+1 << " particle: " << history[i].p() << " quality: " << std::endl;
			//for (unsigned j=0; j<nJets; ++j) {
			//	std::cout << evt.copySequence.clustHist()[i].p().size() << " ";
			//}
			std::cout << std::endl;
		}

		return true;
	}

private:
	typedef typename Event::clust_seq_type::cluster_type cluster_type;

};

#endif // JETHISTORYCOPY_HH_
