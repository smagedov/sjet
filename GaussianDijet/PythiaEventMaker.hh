#ifndef MINIMALEVENTMAKER_HH_
#define MINIMALEVENTMAKER_HH_

#include <random>
#include <cassert>

#include "AbsFrameworkModule.hh"

template <class Event>
class PythiaEventMaker : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline PythiaEventMaker(const std::string& i_label,
                             const std::string& pythiacard)
        : Base(i_label),pythiacard_(pythiacard),firstcall_(true){}
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

        // Return allowing other modules to proceed
        return true;
    }



private:
    //Rng& gen_;
    //std::uniform_real_distribution<> uni_;
    std::shared_ptr<Pythia8::Pythia> pythia_;
    std::string pythiacard_;
    bool firstcall_;
};

#endif // MINIMALEVENTMAKER_HH_
