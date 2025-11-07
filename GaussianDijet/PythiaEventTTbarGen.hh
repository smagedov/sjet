#ifndef PYTHIAEVENTTTBARGEN_HH_
#define PYTHIAEVENTTTBARGEN_HH_

#include <random>
#include <cassert>
#include "Pythia8/Pythia.h"
#include "AbsFrameworkModule.hh"
using namespace Pythia8;
std::vector<Particle> getPartonsFromTTbar(Event &event){
        std::vector<Particle> prompts;
        int iTop = 0;
        int iAntiTop = 0;
        for (int i = 0; i < event.size(); ++i) {
                auto &p = event[i];
                if (p.id() == 6) iTop = i;
                if (p.id() == -6) iAntiTop = i;
        }
        if( event[event[iTop].daughter1()].id() == 5) {
                prompts.push_back(event[event[iTop].daughter1()]);
                int id2 = event[iTop].daughter2();
                while(abs(event[event[id2].daughter1()].id())>6) id2 = event[id2].daughter1();
                prompts.push_back(event[event[id2].daughter1()]);
                prompts.push_back(event[event[id2].daughter2()]);
        }
        else {
                prompts.push_back(event[event[iTop].daughter2()]);
                int id1 = event[iTop].daughter1();
                while(abs(event[event[id1].daughter1()].id())>6) id1 = event[id1].daughter1();
                prompts.push_back(event[event[id1].daughter1()]);
                prompts.push_back(event[event[id1].daughter2()]);

        }
        if( event[event[iAntiTop].daughter1()].id() == -5) {
                prompts.push_back(event[event[iAntiTop].daughter1()]);
                int id2 = event[iAntiTop].daughter2();
                while(abs(event[event[id2].daughter1()].id())>6) id2 = event[id2].daughter1();
                prompts.push_back(event[event[id2].daughter1()]);
                prompts.push_back(event[event[id2].daughter2()]);
        }
        else {
                prompts.push_back(event[event[iAntiTop].daughter2()]);
                int id1 = event[iAntiTop].daughter1();
                while(abs(event[event[id1].daughter1()].id())>6) id1 = event[id1].daughter1();
                prompts.push_back(event[event[id1].daughter1()]);
                prompts.push_back(event[event[id1].daughter2()]);

        }
        return prompts;

}

std::vector<int> getPartonIndicesFromTTbar(Event &event){
        std::vector<int> indices;
        int iTop = 0;
        int iAntiTop = 0;
        for (int i = 0; i < event.size(); ++i) {
                auto &p = event[i];
                if (p.id() == 6) iTop = i;
                if (p.id() == -6) iAntiTop = i;
        }
        if( event[event[iTop].daughter1()].id() == 5) {
                indices.push_back(event[iTop].daughter1());
                int id2 = event[iTop].daughter2();
                while(abs(event[event[id2].daughter1()].id())>6) id2 = event[id2].daughter1();
                indices.push_back(event[id2].daughter1());
                indices.push_back(event[id2].daughter2());
        }
        else {
                indices.push_back(event[iTop].daughter2());
                int id1 = event[iTop].daughter1();
                while(abs(event[event[id1].daughter1()].id())>6) id1 = event[id1].daughter1();
                indices.push_back(event[id1].daughter1());
                indices.push_back(event[id1].daughter2());

        }
        if( event[event[iAntiTop].daughter1()].id() == -5) {
                indices.push_back(event[iAntiTop].daughter1());
                int id2 = event[iAntiTop].daughter2();
                while(abs(event[event[id2].daughter1()].id())>6) id2 = event[id2].daughter1();
                indices.push_back(event[id2].daughter1());
                indices.push_back(event[id2].daughter2());
        }
        else {
                indices.push_back(event[iAntiTop].daughter2());
                int id1 = event[iAntiTop].daughter1();
                while(abs(event[event[id1].daughter1()].id())>6) id1 = event[id1].daughter1();
                indices.push_back(event[id1].daughter1());
                indices.push_back(event[id1].daughter2());

        }
        return indices;

}
void getStableDescendantsRecursive(Event& event, int iParton,std::vector<int>& products) {
    const auto& p = event[iParton];

    for (int iDau = p.daughter1(); iDau <= p.daughter2(); ++iDau) {
        if (iDau <= 0 || iDau >= event.size()) continue;

        const auto& d = event[iDau];

        if (d.isFinal()) {
            products.push_back(iDau);
        }
        else {
            getStableDescendantsRecursive(event, iDau, products);
        }
    }
}

void getStableDescendants(Event& event, int iParton,std::vector<int>& products) {
    const auto& p = event[iParton];
    getStableDescendantsRecursive(event, iParton, products);
    std::sort(products.begin(),products.end());
    auto it_p = std::unique(products.begin(),products.end());
    products.erase(it_p, products.end());
}

void getFinalPartonDescendantsRecursive(Event& event, int iParton,std::vector<int>& products) {
    const auto& p = event[iParton];

    for (int iDau = p.daughter1(); iDau <= p.daughter2(); ++iDau) {
        if (iDau <= 0 || iDau >= event.size()) continue;

        const auto& d = event[iDau];

        if (d.isFinalPartonLevel()) {
            products.push_back(iDau);
        }
        else {
            getFinalPartonDescendantsRecursive(event, iDau, products);
        }
    }
}

void getFinalPartonDescendants(Event& event, int iParton,std::vector<int>& products) {
    const auto& p = event[iParton];
    getFinalPartonDescendantsRecursive(event, iParton, products);
    std::sort(products.begin(),products.end());
    auto it_p = std::unique(products.begin(),products.end());
    products.erase(it_p, products.end());
}

bool have_common_elements(std::vector<int>& v1, std::vector<int>& v2) {
        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());

        std::vector<int> intersection;
        std::set_intersection(v1.begin(), v1.end(),
                              v2.begin(), v2.end(),
                              std::back_inserter(intersection));

        return !intersection.empty();
    }

void getJets(Event& event,vector<int> partons, vector<vector<int>>& jets){
        std::vector<int> all_products;
        for(int i=0; i< partons.size();i++){
                const auto& p = event[partons[i]];
                std::vector<int> stableproducts_fp;
                getStableDescendants(event,partons[i],stableproducts_fp);
                if(have_common_elements(all_products,stableproducts_fp)){
                        continue;
                }
                else{
                        all_products.insert(all_products.end(),stableproducts_fp.begin(),stableproducts_fp.end());
                        jets.push_back(stableproducts_fp);
                        std::cout<<"jet ";
                        for(int j=0;j<stableproducts_fp.size();j++){
                                cout<<stableproducts_fp[j]<<" ";
                        }
                        std::cout<<std::endl;
                }
        }
}

template <class Event>
class PythiaEventTTbarGen : public frw::AbsFrameworkModule<Event>
{
public:
    typedef Event event_type;
    typedef frw::AbsFrameworkModule<Event> Base;

    inline PythiaEventTTbarGen(const std::string& i_label)
        : Base(i_label){}
    inline virtual ~PythiaEventTTbarGen() override {}

    inline virtual PythiaEventTTbarGen* clone() const override
        {return new PythiaEventTTbarGen(*this);}

    inline virtual bool process(Event& evt) override
    {   
        // Make sure that the event has been initialized
        assert(evt.pythiaEventReady);
	assert(!evt.genJetsReady);
	std::vector<int> prompt_indices = getPartonIndicesFromTTbar(*evt.pythiaEvent);
	for(int i=0; i < 6; i++){
		std::vector<int> finalpartons;
		getFinalPartonDescendants(*evt.pythiaEvent, prompt_indices[i], finalpartons);
		evt.promptPartons.insert(evt.promptPartons.end(),finalpartons.begin(),finalpartons.end());
	}
	getJets(*evt.pythiaEvent,evt.promptPartons, evt.genClusters);
	//cout<<"debug "<< "prompt_indices.size()"<<prompt_indices.size()<<endl;
	//cout<<"debug evt.pythiaEvent->size() "<<evt.pythiaEvent->size()<<endl;
        evt.genJetsReady = true;
        // Return allowing other modules to proceed
        return true;
    }

};

#endif // PYTHIAEVENTTTBARGEN_HH_
