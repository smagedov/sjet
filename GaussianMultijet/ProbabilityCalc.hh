#ifndef PROBABILITYCALC_HH_
#define PROBABILITYCALC_HH_

#include <cmath>

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

		const std::vector<cluster_type>& history = evt.diffusionSequence.clustHist();
		const int size = history.size();

		for (int i=0; i<size; ++i) {
			evt.diffusionSequence.setProbability(i, 0.5);
		}

		for (int i=0; i<size; ++i) {
			std::cout << evt.diffusionSequence.getLogProbability(i) << std::endl;
		}

		return true;
	}

private:
	typedef typename Event::clust_seq_type::cluster_type cluster_type;
};

#endif // PROBABILITYCALC_HH_
