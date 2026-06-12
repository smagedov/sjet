#ifndef JETHISTORYCOPY_HH_
#define JETHISTORYCOPY_HH_

#include <cmath>
#include <string>

#include "cnpy.h"
#include "ParticleDistances.hh"


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
		//const unsigned nJets = evt.genJets.size();

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

		std::vector<std::vector<rk::P4>> histcopy(history.size(), std::vector<rk::P4>(nGenClus));

		for (unsigned i=0; i<nParts[nGenClus-1]; ++i) {
			for (unsigned j=0; j<nGenClus; ++j) {
				if (i < nParts[j]) {
					histcopy[i][j] = history[i].p();
					break;
				}
			}
		}

		for (unsigned i=nParts[nGenClus-1]+1; i<size; ++i) {
			int pind1 = history[i].parent1();
			int pind2 = history[i].parent2();
			for (unsigned j=0; j<nGenClus; ++j) {
				histcopy[i][j] = histcopy[pind1][j] + histcopy[pind2][j];
			}
		}

		std::vector<double> leadingPt;
		std::vector<unsigned> leadingId;
		for (unsigned i = 0; i<size; ++i) {
			std::cout << history[i].prob() << std::endl;
			double leadPt = 0.0;
			unsigned leadjetid = 0;
			for (unsigned j = 0; j<nGenClus; ++j) {
				double tmp = histcopy[i][j].pt();
				if (tmp > leadPt) {
					leadPt = tmp;
					leadjetid = j;
				}
			}
			leadingPt.push_back(leadPt);
			leadingId.push_back(leadjetid);
		}

		std::vector<unsigned> labels;
		for (unsigned i = nParts[nGenClus-1]; i<size-1; ++i) {
			unsigned label = 0;
			int pind1 = history[i].parent1();
			int pind2 = history[i].parent2();
			if ((leadingId[pind1] == leadingId[i]) and (leadingId[pind2] == leadingId[i])) {
				label = 1;
			}
			labels.push_back(label);
		}

		const ParticleDeltaR dRcalculator;

		std::vector<double> dists;
		std::vector<double> ratios;
		std::vector<double> massp1;
		std::vector<double> massp2;
		std::vector<double> massd;
		std::vector<double> pt1;
		std::vector<double> eta1;
		std::vector<double> phi1;
		std::vector<double> pt2;
		std::vector<double> eta2;
		std::vector<double> phi2;
		std::vector<double> ptd;
		std::vector<double> deltar;
		for (unsigned i=nParts[nGenClus-1]; i<size-1; ++i) {
			int pind1 = history[i].parent1();
			int pind2 = history[i].parent2();
			int dind = history[i].daughter();
			double dR = dRcalculator(history[pind1].p(), history[pind2].p());
			deltar.push_back(dR);
			ratios.push_back(leadingPt[i]/history[i].scalarPtSum());
			dists.push_back(history[i].dist());
			massp1.push_back(history[pind1].p().m());
			massp2.push_back(history[pind2].p().m());
			massd.push_back(history[dind].p().m());
			pt1.push_back(history[pind1].p().pt());
			eta1.push_back(history[pind1].p().eta());
			phi1.push_back(history[pind1].p().phi());
			pt2.push_back(history[pind2].p().pt());
                        eta2.push_back(history[pind2].p().eta());
                        phi2.push_back(history[pind2].p().phi());
			ptd.push_back(history[dind].p().pt());

		}

		std::string evtnum = std::to_string(evt.number());

		cnpy::npy_save("npyarrays/deltar/" + evtnum + "_deltar.npy", &deltar[0], {deltar.size()}, "w");
		cnpy::npy_save("npyarrays/ratios/" + evtnum + "_ratios.npy", &ratios[0], {ratios.size()}, "w");
		cnpy::npy_save("npyarrays/dists/" + evtnum + "_dists.npy", &dists[0], {dists.size()}, "w");
		cnpy::npy_save("npyarrays/massp1/" + evtnum + "_massp1.npy", &massp1[0], {massp1.size()}, "w");
		cnpy::npy_save("npyarrays/massp2/" + evtnum + "_massp2.npy", &massp2[0], {massp2.size()}, "w");
		cnpy::npy_save("npyarrays/massd/" + evtnum + "_massd.npy", &massd[0], {massd.size()}, "w");
		cnpy::npy_save("npyarrays/pt1/" + evtnum + "_pt1.npy", &pt1[0], {pt1.size()}, "w");
		cnpy::npy_save("npyarrays/eta1/" + evtnum + "_eta1.npy", &eta1[0], {eta1.size()}, "w");
		cnpy::npy_save("npyarrays/phi1/" + evtnum + "_phi1.npy", &phi1[0], {phi1.size()}, "w");
		cnpy::npy_save("npyarrays/pt2/" + evtnum + "_pt2.npy", &pt2[0], {pt2.size()}, "w");
		cnpy::npy_save("npyarrays/eta2/" + evtnum + "_eta2.npy", &eta2[0], {eta2.size()}, "w");
		cnpy::npy_save("npyarrays/phi2/" + evtnum + "_phi2.npy", &phi2[0], {phi2.size()}, "w");
		cnpy::npy_save("npyarrays/ptd/" + evtnum + "_ptd.npy", &ptd[0], {ptd.size()}, "w");
		cnpy::npy_save("npyarrays/labels/" + evtnum + "_labels.npy", &labels[0], {labels.size()}, "w");

		//for (unsigned i = 0; i<history.size(); ++i) {
			//auto& vec = evt.copySequence.clustHist()[i].p();
			//vec.insert(vec.end(), histcopy[i].begin(), histcopy[i].end());
			//std::cout << "evtHist: " <<  evt.diffusionSequence.clustHist()[i].p() << " copyHist: ";
			//std::cout << "i: " << i+1 << " particle: " << history[i].p() << " quality: " << std::endl;
			//for (unsigned j=0; j<nJets; ++j) {
			//	std::cout << evt.copySequence.clustHist()[i].p().size() << " ";
			//}
			//std::cout << std::endl;
		//}
		
		//std::vector<double> jetr;
		//std::vector<double> jetdis;
		//double alpha = 1.5;
		//for (unsigned i=nParts[nGenClus-1]; i<size; ++i) {
		//	unsigned closesti = 0;
		//	double closestdist = dRcalculator(evt.genJets[0], history[i].p());
		//	for (unsigned j=1; j<nJets; ++j) {
		//		double tmp = dRcalculator(evt.genJets[j], history[i].p());
		//		if (tmp < closestdist) {
		//			closestdist = tmp;
		//			closesti = j;
		//		}
		//	}
		//	double tmpclusdis = closestdist + alpha*log(abs(evt.genJets[closesti].pt()/history[i].p().pt()));
		//	jetr.push_back(closestdist);
		//	jetdis.push_back(tmpclusdis);
		//}
		//cnpy::npy_save("npyarrays/jetr/" + evtnum + "_jetr.npy", &jetr[0], {jetr.size()}, "w");
		//cnpy::npy_save("npyarrays/jetdis/" + evtnum + "_jetdis.npy", &jetdis[0], {jetdis.size()}, "w");

		//std::vector<unsigned> mask(size);
		//std::vector<unsigned> closestjeti;
		//for (unsigned i=0; i<nJets; ++i) {
		//	unsigned closesti = 0;
		//	double closestdist = dRcalculator(evt.genJets[i], history[nParts[nGenClus-1]].p()) + alpha*log(abs(evt.genJets[i].pt()/history[nParts[nGenClus-1]].p().pt()));
		//	for (unsigned j=nParts[nJets]; j<size; ++j) {
		//		if (mask[j] == 0) {
		//			double tmp = dRcalculator(evt.genJets[i], history[j].p()) + + alpha*log(abs(evt.genJets[i].pt()/history[j].p().pt()));
		//			if (tmp < closestdist) {
		//				closestdist = tmp;
		//				closesti = j;
		//			}
		//		}
		//	}
		//	closestjeti.push_back(closesti);
		//	mask[closesti] = 1;
		//	int parent1 = history[closesti].parent1();
		//	int parent2 = history[closesti].parent2();
		//	if (parent1 != -1) {
		//		mask[parent1] = 1;
		//	}
		//	if (parent2 != -1) {
		//		mask[parent2] = 1;
		//	}
		//	std::cout << "parent1 id: " << history[closesti].parent1() << " parent2 id: " << history[closesti].parent2() << std::endl;
		//}

		dists.clear();
		ratios.clear();
		massp1.clear();
		massp2.clear();
		massd.clear();
		pt1.clear();
		eta1.clear();
		phi1.clear();
		pt2.clear();
		eta2.clear();
		phi2.clear();
		ptd.clear();
		deltar.clear();
		//jetr.clear();
		//jetdis.clear();
		//mask.clear();
		//closestjeti.clear();

                std::cout << "HOLA" << std::endl;

		return true;
	}

private:
	typedef typename Event::clust_seq_type::cluster_type cluster_type;

};

#endif // JETHISTORYCOPY_HH_
