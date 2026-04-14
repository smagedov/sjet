#ifndef PYTHIAEVENTMAKER_HH_
#define PYTHIAEVENTMAKER_HH_

#include <random>
#include <cassert>

#include "AbsFrameworkModule.hh"
#include "MultiParticleCollectionMaker.hh"
#include "filterByPt.hh"

int poisson(double nAvg, Pythia8::Rndm& rndm) {

  // Set maximum to avoid overflow.
  const int NMAX = 100;

  // Random number.
  double rPoisson = rndm.flat() * exp(nAvg);

  // Initialize.
  double rSum  = 0.;
  double rTerm = 1.;

  // Add to sum and check whether done.
  for (int i = 0; i < NMAX; ) {
    rSum += rTerm;
    if (rSum > rPoisson) return i;

    // Evaluate next term.
    ++i;
    rTerm *= nAvg / i;
  }

  // Emergency return.
  return NMAX;
}

template <class Event>
class PythiaEventMaker : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline PythiaEventMaker(const std::string& i_label,
                             const std::string& pythiacard)
        //With Pileup
	//: Base(i_label),pythiacard_(pythiacard),pythiaPUcard_("MB.cmnd"),firstcall_(true){}
	//Without Pileup
	: Base(i_label),pythiacard_(pythiacard),firstcall_(true){}
    void init(){
	    pythia_ = std::shared_ptr<Pythia8::Pythia>(new Pythia8::Pythia());

	    pythia_->readFile(pythiacard_);
	pythia_->init();
	//pythiaPU_ = std::shared_ptr<Pythia8::Pythia>(new Pythia8::Pythia());
	//    pythiaPU_->readFile(pythiaPUcard_);
	//    pythiaPU_->init();
    }
    inline virtual ~PythiaEventMaker() override {}

    inline virtual PythiaEventMaker* clone() const override
        {return new PythiaEventMaker(*this);}

    inline virtual bool process(Event& evt) override
    {   
	if(firstcall_){
		init();
		firstcall_=false;
	}
        // Make sure that the event has been initialized
        assert(evt.isNumberValid());

	//double nPileupAvg = 1;
        // Make sure that the event has been initialized
	std::cout<<"next "<<pythia_->next()<<std::endl;
	evt.hardScatterSize = pythia_->event.size();
	//int nPileup = poisson(nPileupAvg, pythiaPU_->rndm);
	//std::cout<<"nPU "<<nPileup<<std::endl;
	//for (int iPileup = 0; iPileup < nPileup; ++iPileup) {
	//	pythiaPU_->next();
	//	std::cout<<"before adding "<<pythia_->event.size()<<std::endl;
	//	std::cout<<"pythiaPU_->size() "<<pythiaPU_->event.size()<<std::endl;
	//	(pythia_->event) += (pythiaPU_->event);
	//	std::cout<<"after adding "<<pythia_->event.size()<<std::endl;
	//}

	// Make sure that we are not overwriting existing data
        assert(!evt.genEventReady);

	evt.pythiaEvent = &pythia_->event;
	std::cout<<"next "<<pythia_->next()<<std::endl;;
	std::cout<<"evt.pythiaEvent->size() "<<evt.pythiaEvent->size()<<std::endl;
        evt.pythiaEventReady = true;

        // Return allowing other modules to proceed
        return true;
    }



private:
    typedef typename Event::particle_type MyParticle;

    //std::uniform_real_distribution<> uni_;
    std::string pythiacard_;
    std::shared_ptr<Pythia8::Pythia> pythia_;
    std::shared_ptr<Pythia8::Pythia> pythiaPU_;
    std::string pythiaPUcard_;
    bool firstcall_;
};

#endif // PYTHIAEVENTMAKER_HH_
