#ifndef PYTHIAEVENTTTBARGEN_HH_
#define PYTHIAEVENTTTBARGEN_HH_

#include <random>
#include <cassert>
#include "rk/rk.hh"
#include "Pythia8/Pythia.h"
#include "ParticleMaker.hh"
#include "AbsFrameworkModule.hh"
using namespace Pythia8;
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

void getAncestorsRecursive(Event& event, int iParticle,std::vector<int>& products) {
    const auto& p = event[iParticle];
    if (p.isFinalPartonLevel()){
            products.push_back(iParticle);
    }
    std::vector<int> mo;
    if(p.mother1() >0 ) mo.push_back(p.mother1());
    if(p.mother2() >0 ) mo.push_back(p.mother2());
    for (int i = 0; i < (int)mo.size(); ++i) {
        int iMo = mo[i];
        if (iMo <= 0 || iMo >= event.size()) continue;
        const auto& m = event[iMo];
        if (m.isFinalPartonLevel()) {
            products.push_back(iMo);
        }
        else {
            getAncestorsRecursive(event, iMo, products);
        }
    }
}

void getAncestors(Event& event, int iParticle,std::vector<int>& products) {
    getAncestorsRecursive(event, iParticle, products);
    std::sort(products.begin(),products.end());
    auto it_p = std::unique(products.begin(),products.end());
    products.erase(it_p, products.end());
}


bool isInList(int index, std::vector<int> list){
	for(long unsigned int i=0; i< list.size(); i++){
		if(index == list[i]) return true;
	}
	return false;
}

int findVectorInList(std::vector<int> indices, std::vector<std::vector<int>> list){
        for(long unsigned int i=0; i< list.size(); i++){
		if(indices.size() == list[i].size()){
			bool find = true;
			for(long unsigned int j=0; j< list[i].size(); j++){
				if(indices[j] != list[i][j]) find = false;
			}
			if(find) return i;
		}
        }
        return -1;
}

void getPromptAncestorsRecursive(Event& event, int iParticle,std::vector<int>& products) {
    const auto& p = event[iParticle];
    if (isInList(iParticle,getPartonIndicesFromTTbar(event))){
            products.push_back(iParticle);
    }
    std::vector<int> mo;
    mo.push_back(p.mother1());
    mo.push_back(p.mother2());
    for (int i = 0; i <= 1; ++i) {
        int iMo = mo[i];
        if (iMo <= 0 || iMo >= event.size()) continue;
        if (isInList(iMo,getPartonIndicesFromTTbar(event))) {
            products.push_back(iMo);
        }
        else {
            getPromptAncestorsRecursive(event, iMo, products);
        }
    }
}

void getPromptAncestors(Event& event, int iParticle,std::vector<int>& products) {
    getPromptAncestorsRecursive(event, iParticle, products);
    std::sort(products.begin(),products.end());
    auto it_p = std::unique(products.begin(),products.end());
    products.erase(it_p, products.end());
}

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

void printFinalPariclesFromTTBar(Event &event){
	std::vector<int> finals;
	cout<<"Final particles from this event: "<<endl;
	for (int i = 0; i < event.size(); ++i) {
		auto &p = event[i];
		if(p.isFinal()){
			finals.push_back(i);
			cout<<i<<" (status "<<p.status()<<") ";
		}
	}
	cout<<endl;
	cout<<"Final particles ids from this event: "<<endl;
        for(long unsigned int i=0; i< finals.size(); ++i) {
                cout<<event[finals[i]].id()<<" ";
        }
        cout<<endl;
	cout<<"number of final particles "<<finals.size()<<endl;

        std::vector<int> finalpartons;
        cout<<"Final partons from this event: "<<endl;
        for (int i = 0; i < event.size(); ++i) {
                auto &p = event[i];
                if(p.isFinalPartonLevel()){
                        finalpartons.push_back(i);
                        cout<<i<<" ";
                }
        }
        cout<<endl;
        cout<<"number of final partons from this event "<<finalpartons.size()<<endl;
	std::vector<int> finals_noancestor;
	for(long unsigned int i=0; i< finals.size(); ++i) {
		std::vector<int> ancestors;
                getAncestors(event,finals[i],ancestors);
		std::vector<int> promptancestors;
                getPromptAncestors(event,finals[i],promptancestors);
		cout<<"final particle "<<finals[i]<<", ancestors ";
		for(long unsigned int k=0; k<ancestors.size(); k++){
			cout<<ancestors[k]<<" ";
		}
		if(promptancestors.size()>0){
		cout<<"prompt ancestors ";
		for(long unsigned int k=0; k<promptancestors.size(); k++){
                        cout<<promptancestors[k]<<" ";
                }
		cout<<endl;
		}
		else{
			cout<<"no prompt ancestors "<<endl;
		}
		bool find=false;
		for(long unsigned int j=0; j< finalpartons.size(); ++j){
			if(event[finals[i]].isAncestor(finalpartons[j])) find=true;
		}
		if(!find) finals_noancestor.push_back(finals[i]);
	}
	cout<<"Final particles with no ancestor in final partons:"<<endl;
	for(long unsigned int i=0; i< finals_noancestor.size(); ++i){
		cout<<finals_noancestor[i]<<"(id "<<event[finals_noancestor[i]].id()<<", status "<<event[finals_noancestor[i]].status()<<", mother1 "<<event[finals_noancestor[i]].mother1()<<", mother2 "<<event[finals_noancestor[i]].mother2()<<") ";
	}
	cout<<endl;
	cout<<"number of final particles with no ancestor in final partons "<<finals_noancestor.size()<<endl;
}


void getFinalParticleClusters(Event& event, vector<vector<int>>& jets, int hardScatterSize = -1){
	std::cout << "hardScatterSize: " << hardScatterSize << std::endl;
	if (hardScatterSize < 0) hardScatterSize = event.size();
	std::vector<int> finals;
	std::cout << "hardScatterSize: " << hardScatterSize << " event size: " << event.size() << std::endl;
	cout<<"Final particles from this event: "<<endl;
        for (int i = 0; i < hardScatterSize; ++i) {
                auto &p = event[i];
                if(p.isFinal()){
                        finals.push_back(i);
			cout<<i<<" ";
                }
        }
	cout<<endl;
	cout<<"number of final particles "<<finals.size()<<endl;

	vector<vector<int>> ancestors_list;
	for(long unsigned int i=0; i< finals.size(); ++i) {
                std::vector<int> ancestors;
                getAncestors(event,finals[i],ancestors);
		int jet_index = findVectorInList(ancestors, ancestors_list);
		if(jet_index < 0){
			ancestors_list.push_back(ancestors);
			vector<int> new_jet;
			new_jet.push_back(finals[i]);
			jets.push_back(new_jet);
		}
		else{
			jets[jet_index].push_back(finals[i]);
		}
	}
	for(long unsigned int i=0; i< jets.size(); i++){
		cout<<"jet "<<i<<": ";
		for(long unsigned int j=0;j < jets[i].size(); j++){
			cout<<jets[i][j]<<" ";
		}
		cout<<", "<<jets[i].size()<<" particles"<<endl;
	}
}


void getFinalParticleClustersFromTTbar(Event& event, vector<vector<int>>& jets, int hardScatterSize = -1){
	if (hardScatterSize < 0) hardScatterSize = event.size();
        std::vector<int> finals;
        cout<<"Final particles from this event: "<<endl;
        for (int i = 0; i < hardScatterSize; ++i) {
                auto &p = event[i];
                if(p.isFinal()){
                        finals.push_back(i);
                        cout<<i<<" ";
                }
        }
        cout<<endl;
        cout<<"number of final particles "<<finals.size()<<endl;

        vector<vector<int>> ancestors_list;
        for(long unsigned int i=0; i< finals.size(); ++i) {
		std::vector<int> promptancestors;
                getPromptAncestors(event,finals[i],promptancestors);
		if(promptancestors.size()<=0) continue;
                int jet_index = findVectorInList(promptancestors, ancestors_list);
                if(jet_index < 0){
                        ancestors_list.push_back(promptancestors);
                        vector<int> new_jet;
                        new_jet.push_back(finals[i]);
                        jets.push_back(new_jet);
                }
                else{
                        jets[jet_index].push_back(finals[i]);
                }
        }
        for(long unsigned int i=0; i< jets.size(); i++){
                cout<<"TTbar jet "<<i<<": ";
                for(long unsigned int j=0;j < jets[i].size(); j++){
                        cout<<jets[i][j]<<" ";
                }
                cout<<", "<<jets[i].size()<<" particles"<<endl;
        }
}


std::vector<int> printFinalPartonsFromWholeEvent(Event &event){
        std::vector<int> finals;
        cout<<"Final partons from this event: "<<endl;
        for (int i = 0; i < event.size(); ++i) {
                auto &p = event[i];
                if(p.isFinalPartonLevel()){
                        finals.push_back(i);
                        cout<<i<<" ";
                }
        }
        cout<<endl;
	cout<<"Final partons status from this event: "<<endl;
        for(long unsigned int i=0; i< finals.size(); ++i) {
                cout<<event[finals[i]].status()<<" ";
        }
        cout<<endl;
        cout<<"number of final partons from this event "<<finals.size()<<endl;
	return finals;
}


void getStableDescendantsRecursive(Event& event, int iParton,std::vector<int>& products) {
    const auto& p = event[iParton];
    if (p.isFinal()){
    //if (p.status()>0){
            products.push_back(iParton);
    }
    std::vector<int> dau;
    dau.push_back(p.daughter1());
    dau.push_back(p.daughter2());
    //cout<<"particle "<<iParton<<", idau1 "<<p.daughter1()<<", idau2 "<<p.daughter2()<<endl;
    for (int i = 0; i <= 1; ++i) {
	int iDau = dau[i];
        if (iDau <= 0 || iDau >= event.size()) continue;
        //cout<<"daughter particle "<<iDau<<" found"<<endl;
        const auto& d = event[iDau];
        if (d.isFinal()) {
	//if (d.status()>0) {
	    //cout<<"daughter particle "<<iDau<<" is final"<<endl;
            products.push_back(iDau);
        }
        else {
	    //cout<<"daughter particle "<<iDau<<" is not final, continue the chain"<<endl;
            getStableDescendantsRecursive(event, iDau, products);
        }
    }
}

void getStableDescendants(Event& event, int iParton,std::vector<int>& products) {
    getStableDescendantsRecursive(event, iParton, products);
    std::sort(products.begin(),products.end());
    auto it_p = std::unique(products.begin(),products.end());
    products.erase(it_p, products.end());
}

void getFinalPartonDescendantsRecursive(Event& event, int iParton,std::vector<int>& products) {
    const auto& p = event[iParton];
    if (p.isFinalPartonLevel()){
	    products.push_back(iParton);
    }
    std::vector<int> dau;
    dau.push_back(p.daughter1());
    dau.push_back(p.daughter2());
    for (int i = 0; i <= 1; ++i) {
	int iDau = dau[i];
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
        for(long unsigned int i=0; i< partons.size();i++){
                const auto& p = event[partons[i]];
                std::vector<int> stableproducts_fp;
                getStableDescendants(event,partons[i],stableproducts_fp);
                if(have_common_elements(all_products,stableproducts_fp)){
			//for(int j=0;j<stableproducts_fp.size();j++){
                        //        cout<<stableproducts_fp[j]<<" ";
                        //}
			//cout<<"overlapping with existing jets "<<endl;
                        continue;
                }
                else{
			cout<<"final parton "<<partons[i]<< ", id "<<p.id()<<endl;
                        all_products.insert(all_products.end(),stableproducts_fp.begin(),stableproducts_fp.end());
                        jets.push_back(stableproducts_fp);
                        std::cout<<"jet ";
                        for(long unsigned int j=0;j<stableproducts_fp.size();j++){
                                cout<<stableproducts_fp[j]<<" ";
                        }
			std::cout<<std::endl;
			std::cout<<"ids ";
			for(long unsigned int j=0;j<stableproducts_fp.size();j++){
                                cout<<event[stableproducts_fp[j]].id()<<" ";
                        }
                        std::cout<<std::endl;
                }
        }
	cout<<"Jet particles status: "<<endl;
        for(long unsigned int i=0; i< all_products.size(); ++i) {
                cout<<event[all_products[i]].status()<<" ";
        }
        cout<<endl;
	cout<<"Jet particles: "<<endl;
        for(long unsigned int i=0; i< all_products.size(); ++i) {
                cout<<all_products[i]<<" ";
        }
        cout<<endl;
	cout<<"total number of particles in these jets "<<all_products.size()<<endl;
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
	getFinalParticleClusters(*evt.pythiaEvent,evt.genClusters,evt.hardScatterSize);
	getFinalParticleClustersFromTTbar(*evt.pythiaEvent,evt.genClustersFromHardCollision,evt.hardScatterSize);
        evt.genJetsReady = true;

	evt.genEvent.first.clear();
        evt.genEvent.second.clear();

        std::vector<rk::P4> finalParticles;
        finalParticles.reserve(evt.pythiaEvent->size());

        // Loop through particles in the event
	for (long unsigned int i=0; i<evt.genClusters.size(); ++i) {
		rk::P4 genJet;
                if (!evt.genClusters[i].empty()) {
			int jetParts = evt.genClusters[i].size();
                        for (int j=0; j<jetParts; ++j) {
				int pid = evt.genClusters[i][j];
				const Pythia8::Particle& p = (*evt.pythiaEvent)[pid];
				if (p.isFinal()) {
					rk::P4 part = rk::P4(p.pT()*geom3::Vector3(cos(p.phi()), sin(p.phi()), sinh(p.eta())), p.m());
					genJet = genJet + part;
					evt.genEvent.first.push_back(part);
					evt.genEvent.second.push_back(jetParts);
				}
                        }
			evt.genJets.push_back(genJet);
		}
	}

        // Return allowing other modules to proceed
        return true;
    }

private:
        typedef typename Event::particle_type MyParticle;
    	typedef std::pair<std::vector<MyParticle>, std::vector<unsigned> > GenEvent;

};

#endif // PYTHIAEVENTTTBARGEN_HH_
