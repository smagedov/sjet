#ifndef PYTHIAEVENTMAKER_HH_
#define PYTHIAEVENTMAKER_HH_

#include <random>
#include <cassert>

#include "AbsFrameworkModule.hh"
#include "MultiParticleCollectionMaker.hh"
#include "filterByPt.hh"

template <class Event, class Rng>
class PythiaEventMaker : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef Rng rng_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline PythiaEventMaker(const std::string& i_label,
		             Rng& gen,
                             const std::string& pythiacard)
        : Base(i_label), gen_(gen), pythiacard_(pythiacard),firstcall_(true){}
    void init(){
	    pythia_ = std::shared_ptr<Pythia8::Pythia>(new Pythia8::Pythia());

	    pythia_->readFile(pythiacard_);
	pythia_->init();
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

	// Make sure that we are not overwriting existing data
        assert(!evt.genEventReady);

	evt.pythiaEvent = &pythia_->event;
	std::cout<<"next "<<pythia_->next()<<std::endl;;
	std::cout<<"evt.pythiaEvent->size() "<<evt.pythiaEvent->size()<<std::endl;
	//while(!evt.pythia->next()) continue;
        //evt.pythiaEvent = &evt.pythia->event;
        evt.pythiaEventReady = true;
        // Make sure that we are not overwriting existing data
        // assert(!evt.randomDatumReady);

        // Fill out this module's portion of event data
        //evt.randomDatum = uni_(gen_);

        // Note that this portion of event is ready
        //evt.randomDatumReady = true;
	const GenEvent& genEvent0 = eventBuilder_.make(gen_);
	
	evt.genEventReady = true;

        // Return allowing other modules to proceed
        return true;
    }



private:
    typedef typename Event::particle_type MyParticle;
    typedef MultiParticleCollectionMaker<MyParticle,Rng> MyEventBuilder;
    typedef std::pair<std::vector<MyParticle>, std::vector<unsigned> > GenEvent;

    //Rng& gen_;
    //std::uniform_real_distribution<> uni_;
    std::shared_ptr<Pythia8::Pythia> pythia_;
    Rng& gen_;
    std::string pythiacard_;
    bool firstcall_;
    MyEventBuilder eventBuilder_;
};

#endif // PYTHIAEVENTMAKER_HH_
