#ifndef PROBABILITYCALC_HH_
#define PROBABILITYCALC_HH_

#include <torch/script.h>
#include <torch/torch.h>
#include <vector>
#include <cmath>

#include "ParticleDistances.hh"

template <class Event>
class ProbabilityCalc : public frw::AbsFrameworkAnalyzer<Event>
{
public:
        typedef Event event_type;
        typedef frw::AbsFrameworkAnalyzer<Event> Base;

        inline ProbabilityCalc(const std::string& i_label)
                : Base(i_label)
        {
        }

        inline virtual ~ProbabilityCalc() override {}

        inline virtual ProbabilityCalc* clone() const override
                {return new ProbabilityCalc(*this);}

        inline virtual bool analyze(const Event& evt) override
        {
		assert(evt.diffusionSequenceReady);
		assert(evt.simpleDiffusionJetsReady);

		try {
			model = torch::jit::load("model.pt");
			model.eval();
		}
		catch (const c10::Error& e) {
			std::cerr << "Error loading TorchScript model\n";
		}

                std::vector<unsigned> nParts;
                unsigned curParts = 0;
                const unsigned nGenClus = evt.genClusters.size();
                for (unsigned i=0; i<nGenClus; ++i) {
                        if (!evt.genClusters[i].empty()) {
                                curParts = curParts + evt.genClusters[i].size();
                                nParts.push_back(curParts);
                        }
                }

		const std::vector<cluster_type>& history = evt.diffusionSequence.clustHist();
		const int size = history.size();

		for (unsigned i=0; i<nParts[nGenClus-1]; ++i) {
			evt.diffusionSequence.setProbability(i, 1);
		}

		const ParticleDeltaR dRcalculator;

		for (unsigned i=nParts[nGenClus-1]; i<size-1; ++i) {
			int pind1 = history[i].parent1();
                        int pind2 = history[i].parent2();
                        int dind = history[i].daughter();
                        double dR = dRcalculator(history[pind1].p(), history[pind2].p());
			std::vector<double> features = {
				history[dind].p().pt(),
			 	history[dind].p().m(),
			 	history[pind1].p().m(),
			 	history[pind2].p().m(),
			 	history[pind1].p().pt(),
			 	history[pind2].p().pt(),
			 	history[i].dist(),
			 	dR
			};

			torch::Tensor input = torch::tensor(features).unsqueeze(0);
			std::vector<torch::jit::IValue> inputs;
                	inputs.push_back(input);

                	torch::Tensor output = model.forward(inputs).toTensor();

                	double prob = output.item<double>();

                	evt.diffusionSequence.setProbability(i, prob);
		}

		//for (int i=0; i<size; ++i) {
		//	evt.diffusionSequence.setProbability(i, 0.5);
		//}

		for (int i=0; i<size; ++i) {
			std::cout << evt.diffusionSequence.getLogProbability(i) << std::endl;
		}

		return true;
	}

private:
	typedef typename Event::clust_seq_type::cluster_type cluster_type;
	torch::jit::script::Module model;
};

#endif // PROBABILITYCALC_HH_
