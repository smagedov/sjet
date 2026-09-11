#ifndef PYTHIAEVENTTTBARGEN_HH_
#define PYTHIAEVENTTTBARGEN_HH_

#include <random>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <set>
#include <limits>
#include "rk/rk.hh"
#include "Pythia8/Pythia.h"
#include "ParticleMaker.hh"
#include "AbsFrameworkModule.hh"
#include "cnpy.h"
using namespace Pythia8;
namespace {

bool validEventIndex(const Event& event, int i) {
    return i > 0 && i < event.size();
}

std::vector<int> daughterIndices(const Event& event, int iParticle) {
    if (!validEventIndex(event, iParticle)) return {};

    std::vector<int> daughters = event[iParticle].daughterList();
    daughters.erase(
        std::remove_if(daughters.begin(), daughters.end(),
                       [&event](int i) { return !validEventIndex(event, i); }),
        daughters.end());
    return daughters;
}

// Follow Pythia copies of the same particle until the copy that actually decays.
int lastCopy(Event& event, int index) {
    if (!validEventIndex(event, index)) return -1;

    std::set<int> visited;
    int current = index;

    while (validEventIndex(event, current) && visited.insert(current).second) {
        int next = -1;
        for (int d : daughterIndices(event, current)) {
            if (event[d].id() == event[current].id()) {
                next = d;
                break;
            }
        }
        if (next < 0) break;
        current = next;
    }
    return current;
}

int findLastParticleWithId(Event& event, int id) {
    int found = -1;
    for (int i = 1; i < event.size(); ++i) {
        if (event[i].id() == id) found = i;
    }
    if (found < 0) return -1;
    return lastCopy(event, found);
}

bool findTopDecay(Event& event, int topIndex,
                  int& bIndex, int& wIndex) {
    bIndex = -1;
    wIndex = -1;

    int top = lastCopy(event, topIndex);
    if (!validEventIndex(event, top)) return false;

    for (int d : daughterIndices(event, top)) {
        const int absId = std::abs(event[d].id());
        // Keep the b quark directly from the top decay. Following its
        // later same-id copies would lose radiation emitted before that copy.
        if (absId == 5)  bIndex = d;
        if (absId == 24) wIndex = lastCopy(event, d);
    }

    return validEventIndex(event, bIndex) && validEventIndex(event, wIndex);
}

bool findHadronicWDecay(Event& event, int wIndex,
                        int& q1Index, int& q2Index) {
    q1Index = -1;
    q2Index = -1;

    int w = lastCopy(event, wIndex);
    if (!validEventIndex(event, w)) return false;

    std::vector<int> quarks;
    for (int d : daughterIndices(event, w)) {
        const int absId = std::abs(event[d].id());
        // Keep the quark directly from the W decay so its complete shower
        // remains underneath this truth ancestor.
        if (absId >= 1 && absId <= 5)
            quarks.push_back(d);
    }

    std::sort(quarks.begin(), quarks.end());
    quarks.erase(std::unique(quarks.begin(), quarks.end()), quarks.end());

    if (quarks.size() != 2) return false;
    q1Index = quarks[0];
    q2Index = quarks[1];
    return true;
}

void getPromptAncestorsRecursiveImpl(Event& event,
                                     int iParticle,
                                     const std::vector<int>& promptPartons,
                                     std::vector<int>& products,
                                     std::set<int>& visited) {
    if (!validEventIndex(event, iParticle)) return;
    if (!visited.insert(iParticle).second) return;

    if (std::find(promptPartons.begin(), promptPartons.end(), iParticle)
        != promptPartons.end()) {
        products.push_back(iParticle);
        return;
    }

    const auto& p = event[iParticle];

    // Use Pythia's motherList()
    const std::vector<int> mothers = p.motherList();
    for (int iMo : mothers) {
        if (!validEventIndex(event, iMo)) continue;
        getPromptAncestorsRecursiveImpl(event, iMo, promptPartons,
                                        products, visited);
    }
}

}

std::vector<int> getPartonIndicesFromTTbar(Event &event) {
    std::vector<int> indices;

    const int iTop = findLastParticleWithId(event, 6);
    const int iAntiTop = findLastParticleWithId(event, -6);

    if (iTop < 0 || iAntiTop < 0) {
        std::cout << "WARNING: could not find both top and anti-top in event"
                  << std::endl;
        return indices;
    }

    int b = -1, wPlus = -1;
    int bbar = -1, wMinus = -1;

    if (!findTopDecay(event, iTop, b, wPlus) ||
        !findTopDecay(event, iAntiTop, bbar, wMinus)) {
        std::cout << "WARNING: could not identify t -> bW and tbar -> bbarW"
                  << std::endl;
        return indices;
    }

    int q1 = -1, q2 = -1;
    int q3 = -1, q4 = -1;

    if (!findHadronicWDecay(event, wPlus, q1, q2) ||
        !findHadronicWDecay(event, wMinus, q3, q4)) {
        std::cout << "WARNING: TTbar event is not fully hadronic, or the W decay "
                     "quarks could not be identified"
                  << std::endl;
        return indices;
    }

    // Fixed ordering:
    //   0 = b from t
    //   1,2 = quarks from W+/- belonging to t
    //   3 = bbar from tbar
    //   4,5 = quarks from the W belonging to tbar
    indices = {b, q1, q2, bbar, q3, q4};
    return indices;
}

void getAncestorsRecursive(Event& event, int iParticle, std::vector<int>& products) {
    if (iParticle <= 0 || iParticle >= event.size()) return;

    const auto& p = event[iParticle];

    // A final parton is the shower endpoint we want to use as the cluster
    // ancestor. Do not continue farther up the history from here.
    if (p.isFinalPartonLevel()) {
        products.push_back(iParticle);
        return;
    }

    const int mothers[2] = {p.mother1(), p.mother2()};
    for (int iMo : mothers) {
        if (iMo <= 0 || iMo >= event.size()) continue;
        getAncestorsRecursive(event, iMo, products);
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

void getPromptAncestorsRecursive(Event& event, int iParticle,
                                 std::vector<int>& products) {
    const std::vector<int> promptPartons = getPartonIndicesFromTTbar(event);
    if (promptPartons.size() != 6) return;

    std::set<int> visited;
    getPromptAncestorsRecursiveImpl(event, iParticle, promptPartons,
                                    products, visited);
}

void getPromptAncestors(Event& event, int iParticle,
                        std::vector<int>& products) {
    const std::vector<int> promptPartons = getPartonIndicesFromTTbar(event);
    if (promptPartons.size() != 6) return;

    std::set<int> visited;
    getPromptAncestorsRecursiveImpl(event, iParticle, promptPartons,
                                    products, visited);
    std::sort(products.begin(), products.end());
    products.erase(std::unique(products.begin(), products.end()), products.end());
}

std::vector<Particle> getPartonsFromTTbar(Event &event) {
    std::vector<Particle> prompts;
    const std::vector<int> indices = getPartonIndicesFromTTbar(event);
    prompts.reserve(indices.size());
    for (int index : indices)
        prompts.push_back(event[index]);
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
    jets.clear();
	std::cout << "hardScatterSize: " << hardScatterSize << std::endl;
	if (hardScatterSize < 0) hardScatterSize = event.size();
	std::vector<int> finals;
	std::cout << "hardScatterSize: " << hardScatterSize << " event size: " << event.size() << std::endl;
	cout<<"Final particles from this event: "<<endl;
        for (int i = 0; i < event.size(); ++i) {
                auto &p = event[i];
                if (p.isFinal() && !p.isFinalPartonLevel()){
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
        if (ancestors.empty()) continue;
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


bool hadronizationSystemMixesTTbarAndNonTTbarImpl(
        Event& event,
        int iParticle,
        const std::vector<int>& ttbarPartons,
        std::set<int>& visited) {
    if (!validEventIndex(event, iParticle)) return false;
    if (!visited.insert(iParticle).second) return false;

    const Particle& p = event[iParticle];
    const int absStatus = std::abs(p.status());

    // Primary hadrons from hadronization.
    if (absStatus >= 81 && absStatus <= 89) {
        const std::vector<int> mothers = p.motherList();

        bool hasTTbarContributor = false;
        bool hasNonTTbarContributor = false;

        for (int iMo : mothers) {
            if (!validEventIndex(event, iMo)) continue;

            // Treat every entry in the expanded motherList() as a contributor
            // to the hadronization system.  Do not pre-filter to q/g only:
            // beam-remnant diquarks and other non-q/g contributors are exactly
            // the objects that can reveal that a hadron is not owned by TTbar.
            std::vector<int> promptAncestors;
            std::set<int> ancestryVisited;
            getPromptAncestorsRecursiveImpl(event, iMo, ttbarPartons,
                                            promptAncestors, ancestryVisited);
            std::sort(promptAncestors.begin(), promptAncestors.end());
            promptAncestors.erase(
                std::unique(promptAncestors.begin(), promptAncestors.end()),
                promptAncestors.end());

            if (promptAncestors.empty())
                hasNonTTbarContributor = true;
            else
                hasTTbarContributor = true;

            if (hasTTbarContributor && hasNonTTbarContributor)
                return true;
        }

        return false;
    }

    for (int iMo : p.motherList()) {
        if (!validEventIndex(event, iMo)) continue;
        if (hadronizationSystemMixesTTbarAndNonTTbarImpl(
                event, iMo, ttbarPartons, visited))
            return true;
    }

    return false;
}

bool hadronizationSystemMixesTTbarAndNonTTbar(
        Event& event,
        int iParticle,
        const std::vector<int>& ttbarPartons) {
    std::set<int> visited;
    return hadronizationSystemMixesTTbarAndNonTTbarImpl(
        event, iParticle, ttbarPartons, visited);
}


bool hasHadronizationAncestorImpl(Event& event, int iParticle, std::set<int>& visited) {
    if (!validEventIndex(event, iParticle)) return false;
    if (!visited.insert(iParticle).second) return false;

    const Particle& p = event[iParticle];
    const int absStatus = std::abs(p.status());
    if (absStatus >= 81 && absStatus <= 89) return true;

    for (int iMo : p.motherList()) {
        if (!validEventIndex(event, iMo)) continue;
        if (hasHadronizationAncestorImpl(event, iMo, visited)) return true;
    }
    return false;
}

bool hasHadronizationAncestor(Event& event, int iParticle) {
    std::set<int> visited;
    return hasHadronizationAncestorImpl(event, iParticle, visited);
}

double deltaPhiTruth(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi > M_PI)  dphi -= 2.0 * M_PI;
    while (dphi <= -M_PI) dphi += 2.0 * M_PI;
    return dphi;
}

double deltaRTruth(const Particle& a, const Particle& b) {
    const double deta = a.eta() - b.eta();
    const double dphi = deltaPhiTruth(a.phi(), b.phi());
    return std::sqrt(deta*deta + dphi*dphi);
}


void collectPrimaryHadronizationAncestorsImpl(
        Event& event,
        int iParticle,
        std::vector<int>& primaryHadrons,
        std::set<int>& visited) {
    if (!validEventIndex(event, iParticle)) return;
    if (!visited.insert(iParticle).second) return;

    const Particle& p = event[iParticle];
    const int absStatus = std::abs(p.status());

    // Status 81--89 entries are primary hadrons produced by hadronization.
    // Stop here: their motherList() describes the partonic hadronization system.
    if (absStatus >= 81 && absStatus <= 89) {
        primaryHadrons.push_back(iParticle);
        return;
    }

    for (int iMo : p.motherList()) {
        if (!validEventIndex(event, iMo)) continue;
        collectPrimaryHadronizationAncestorsImpl(
            event, iMo, primaryHadrons, visited);
    }
}

std::vector<int> collectPrimaryHadronizationAncestors(
        Event& event, int iParticle) {
    std::vector<int> primaryHadrons;
    std::set<int> visited;
    collectPrimaryHadronizationAncestorsImpl(
        event, iParticle, primaryHadrons, visited);
    std::sort(primaryHadrons.begin(), primaryHadrons.end());
    primaryHadrons.erase(
        std::unique(primaryHadrons.begin(), primaryHadrons.end()),
        primaryHadrons.end());
    return primaryHadrons;
}

// For a stable particle that ultimately came from hadronization, identify the
// geometrically closest contributor to the primary hadronization system.
// Pythia's motherList() expands mother ranges, so a primary hadron with
// mothers (i1,i2) can expose all contributors in that range.
//
// IMPORTANT: do not pre-filter the contributors by PDG id.  In particular,
// beam-remnant diquarks or other non-q/g entries must be allowed to win the
// geometric comparison.  Only after the closest contributor is chosen do we
// ask whether that contributor can be traced back to one of the six TTbar
// decay quarks.
int closestHadronizationParent(Event& event,
                               int iParticle,
                               const std::vector<int>& primaryHadrons,
                               double* bestDeltaR = nullptr) {
    if (!validEventIndex(event, iParticle)) return -1;

    const Particle& stable = event[iParticle];
    int bestParent = -1;
    double bestDR = std::numeric_limits<double>::infinity();

    for (int iHadron : primaryHadrons) {
        if (!validEventIndex(event, iHadron)) continue;

        for (int iMo : event[iHadron].motherList()) {
            if (!validEventIndex(event, iMo)) continue;

            // Every actual contributor competes geometrically.  The selected
            // contributor's ancestry is tested only after this choice.
            const double dR = deltaRTruth(stable, event[iMo]);
            if (dR < bestDR) {
                bestDR = dR;
                bestParent = iMo;
            }
        }
    }

    if (bestDeltaR) *bestDeltaR = bestDR;
    return bestParent;
}

// Map one partonic parent back to the six hard TTbar decay quarks.  If the
// selected parent itself has more than one TTbar ancestor, retain the existing
// ambiguity rule and choose the closest hard TTbar ancestor in DeltaR.
int ttbarJetFromSelectedParent(Event& event,
                               int iStableParticle,
                               int iParent,
                               const std::vector<int>& ttbarPartons,
                               bool* wasAmbiguous = nullptr) {
    if (wasAmbiguous) *wasAmbiguous = false;
    if (!validEventIndex(event, iStableParticle) ||
        !validEventIndex(event, iParent))
        return -1;

    std::vector<int> promptAncestors;
    std::set<int> visited;
    getPromptAncestorsRecursiveImpl(event, iParent, ttbarPartons,
                                    promptAncestors, visited);
    std::sort(promptAncestors.begin(), promptAncestors.end());
    promptAncestors.erase(
        std::unique(promptAncestors.begin(), promptAncestors.end()),
        promptAncestors.end());

    if (promptAncestors.empty()) return -1;

    if (promptAncestors.size() == 1) {
        auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(),
                            promptAncestors[0]);
        if (it == ttbarPartons.end()) return -1;
        return static_cast<int>(std::distance(ttbarPartons.begin(), it));
    }

    if (wasAmbiguous) *wasAmbiguous = true;

    const Particle& stable = event[iStableParticle];
    int bestJet = -1;
    double bestDR = std::numeric_limits<double>::infinity();

    for (int ancestor : promptAncestors) {
        auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(), ancestor);
        if (it == ttbarPartons.end()) continue;

        const double dR = deltaRTruth(stable, event[ancestor]);
        if (dR < bestDR) {
            bestDR = dR;
            bestJet = static_cast<int>(
                std::distance(ttbarPartons.begin(), it));
        }
    }

    return bestJet;
}

void getFinalParticleClustersFromTTbar(Event& event,
                                        vector<vector<int>>& jets,
                                        int hardScatterSize = -1,
                                        std::vector<std::size_t>* rejectedMixedPerJet = nullptr) {
    (void)hardScatterSize;
    jets.clear();

    const std::vector<int> ttbarPartons = getPartonIndicesFromTTbar(event);
    if (ttbarPartons.size() != 6) {
        std::cout << "WARNING: expected 6 TTbar decay quarks, found "
                  << ttbarPartons.size() << std::endl;
        return;
    }

    // Fixed six-entry layout, one truth collection per TTbar decay quark.
    jets.resize(6);
    if (rejectedMixedPerJet) rejectedMixedPerJet->assign(6, 0);

    std::size_t ambiguousParticles = 0;
    std::size_t mixedHadronizationParticles = 0;
    std::size_t hadronizationParentAssigned = 0;
    std::size_t hadronizationParentRejectedNonTTbar = 0;
    std::size_t hadronizationParentRejectedFar = 0;
    std::size_t hadronizationParentMissing = 0;

    // This cut is relative to the selected local hadronization contributor,
    // NOT relative to the original hard TTbar quark.
    constexpr double maxHadronizationParentDeltaR = 1.0;

    for (int i = 1; i < event.size(); ++i) {
        const auto& p = event[i];

        // Stable hadron-level particles only; do not admit shower partons.
        if (!p.isFinal() || p.isFinalPartonLevel()) continue;

        const std::vector<int> primaryHadrons =
            collectPrimaryHadronizationAncestors(event, i);

        // Hadronization-origin particles use parent matching rather than a
        // fixed DeltaR cut to the hard TTbar quark.  Geometry chooses the
        // closest contributing parton; ancestry of that selected parton then
        // decides whether the stable particle belongs to one of the six TTbar
        // truth jets.
        if (!primaryHadrons.empty()) {
            const bool mixedTTbarNonTTbar =
                hadronizationSystemMixesTTbarAndNonTTbar(
                    event, i, ttbarPartons);

            if (mixedTTbarNonTTbar) {
                ++mixedHadronizationParticles;
            }

            double parentDeltaR = std::numeric_limits<double>::infinity();
            const int closestParent =
                closestHadronizationParent(event, i, primaryHadrons,
                                           &parentDeltaR);

            if (closestParent < 0) {
                ++hadronizationParentMissing;
                continue;
            }

            bool parentAmbiguous = false;
            const int jetIndex =
                ttbarJetFromSelectedParent(event, i, closestParent,
                                           ttbarPartons, &parentAmbiguous);

            if (parentAmbiguous) ++ambiguousParticles;

            if (jetIndex < 0) {
                ++hadronizationParentRejectedNonTTbar;
                continue;
            }

            // The closest contributor must also be genuinely local to the
            // produced hadron.  This removes string fragments that happen to
            // be closer to the TTbar end of a long string than to the beam
            // remnant, while still allowing wide-angle shower radiation:
            // the comparison is to the local shower/hadronization parent,
            // not to the original hard quark.
            if (parentDeltaR > maxHadronizationParentDeltaR) {
                ++hadronizationParentRejectedFar;
                continue;
            }

            jets[jetIndex].push_back(i);
            ++hadronizationParentAssigned;

            // Preserve the old per-jet mixed-system diagnostic.  It is now
            // purely informational and does not control assignment.
            if (mixedTTbarNonTTbar && rejectedMixedPerJet &&
                jetIndex < static_cast<int>(rejectedMixedPerJet->size())) {
                ++(*rejectedMixedPerJet)[jetIndex];
            }

            continue;
        }

        // No hadronization ancestor: retain the original ancestry rule for
        // perturbative/direct stable descendants.  Multiple TTbar ancestors
        // are resolved geometrically exactly as before.
        std::vector<int> promptAncestors;
        std::set<int> visited;
        getPromptAncestorsRecursiveImpl(event, i, ttbarPartons,
                                        promptAncestors, visited);
        std::sort(promptAncestors.begin(), promptAncestors.end());
        promptAncestors.erase(
            std::unique(promptAncestors.begin(), promptAncestors.end()),
            promptAncestors.end());

        if (promptAncestors.empty()) continue;

        if (promptAncestors.size() == 1) {
            auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(),
                                promptAncestors[0]);
            if (it != ttbarPartons.end()) {
                const std::size_t jetIndex =
                    std::distance(ttbarPartons.begin(), it);
                jets[jetIndex].push_back(i);
            }
            continue;
        }

        ++ambiguousParticles;
        int bestJet = -1;
        double bestDeltaR = std::numeric_limits<double>::infinity();

        for (int ancestor : promptAncestors) {
            auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(), ancestor);
            if (it == ttbarPartons.end()) continue;

            const int jetIndex =
                static_cast<int>(std::distance(ttbarPartons.begin(), it));
            const double dR = deltaRTruth(p, event[ancestor]);

            if (dR < bestDeltaR) {
                bestDeltaR = dR;
                bestJet = jetIndex;
            }
        }

        if (bestJet >= 0)
            jets[bestJet].push_back(i);
    }

    std::cout << "TTbar truth jets (hadronization parent matching; ancestry otherwise):"
              << std::endl;
    for (std::size_t i = 0; i < jets.size(); ++i) {
        const int partonIndex = ttbarPartons[i];
        std::cout << "  TTbar jet " << i
                  << " -> parton " << partonIndex
                  << " (id " << event[partonIndex].id() << ")"
                  << ", " << jets[i].size() << " stable descendants"
                  << std::endl;
    }
    std::cout << "  ambiguous selected-parent/TTbar assignments resolved with DeltaR: "
              << ambiguousParticles << std::endl;
    std::cout << "  particles from mixed TTbar/non-TTbar hadronization systems seen: "
              << mixedHadronizationParticles << std::endl;
    std::cout << "  hadronization-origin particles assigned through closest parent: "
              << hadronizationParentAssigned << std::endl;
    std::cout << "  hadronization-origin particles rejected because closest parent has no TTbar ancestry: "
              << hadronizationParentRejectedNonTTbar << std::endl;
    std::cout << "  hadronization-origin particles rejected because DeltaR(stable, closest parent) > "
              << maxHadronizationParentDeltaR << ": "
              << hadronizationParentRejectedFar << std::endl;
    std::cout << "  hadronization-origin particles with no usable hadronization contributor: "
              << hadronizationParentMissing << std::endl;
}




void printBCollectionSummary(Event& event,
                             const std::vector<int>& collection,
                             int hardBIndex,
                             const char* label,
                             double reconstructedTopMass,
                             std::size_t rejectedMixedCount) {
    if (!validEventIndex(event, hardBIndex)) return;

    double sumPx = 0.0, sumPy = 0.0, sumPz = 0.0, sumE = 0.0;
    double maxAbsEta = 0.0;
    double maxDeltaR = 0.0;
    std::size_t nFinal = 0;

    for (int idx : collection) {
        if (!validEventIndex(event, idx)) continue;
        const auto& p = event[idx];
        if (!p.isFinal()) continue;

        ++nFinal;
        sumPx += p.px();
        sumPy += p.py();
        sumPz += p.pz();
        sumE  += p.e();
        maxAbsEta = std::max(maxAbsEta, std::abs(p.eta()));
        maxDeltaR = std::max(maxDeltaR, deltaRTruth(p, event[hardBIndex]));
    }

    const double pt = std::sqrt(sumPx*sumPx + sumPy*sumPy);
    const double m2 = sumE*sumE - sumPx*sumPx - sumPy*sumPy - sumPz*sumPz;
    const double mass = std::sqrt(std::max(0.0, m2));

    std::cout << "TTBAR_B_SUMMARY " << label
              << " recTopMass " << reconstructedTopMass
              << " nConstituents " << nFinal
              << " pt " << pt
              << " mass " << mass
              << " energy " << sumE
              << " maxAbsEta " << maxAbsEta
              << " maxDeltaR " << maxDeltaR
              << " mixedHadronization " << rejectedMixedCount
              << std::endl;
}

std::vector<int> shortestMotherPathToTarget(Event& event,
                                             int startIndex,
                                             int targetIndex) {
    std::vector<int> empty;
    if (!validEventIndex(event, startIndex) || !validEventIndex(event, targetIndex))
        return empty;

    std::vector<int> previous(event.size(), -2);
    std::vector<int> current{startIndex};
    previous[startIndex] = -1;

    while (!current.empty()) {
        std::vector<int> next;

        for (int idx : current) {
            if (idx == targetIndex) {
                std::vector<int> path;
                int p = idx;
                while (p >= 0) {
                    path.push_back(p);
                    p = previous[p];
                }
                std::reverse(path.begin(), path.end());
                return path;
            }

            const int mothers[2] = {event[idx].mother1(), event[idx].mother2()};
            for (int m : mothers) {
                if (!validEventIndex(event, m)) continue;
                if (previous[m] != -2) continue;
                previous[m] = idx;
                next.push_back(m);
            }
        }

        current.swap(next);
    }

    return empty;
}


void printHadronizationParentDecisionDiagnostic(
        Event& event,
        int stableIndex,
        const std::vector<int>& ttbarPartons) {
    if (!validEventIndex(event, stableIndex)) return;

    const Particle& stable = event[stableIndex];
    const std::vector<int> primaryHadrons =
        collectPrimaryHadronizationAncestors(event, stableIndex);

    if (primaryHadrons.empty()) {
        std::cout << "      closest-parent diagnostic: no primary hadronization ancestor"
                  << std::endl;
        return;
    }

    std::cout << "      closest-parent diagnostic:" << std::endl;

    for (int iHadron : primaryHadrons) {
        if (!validEventIndex(event, iHadron)) continue;

        const Particle& hadron = event[iHadron];
        std::cout << "        primary hadron idx " << iHadron
                  << " id " << hadron.id()
                  << " status " << hadron.status()
                  << " eta " << hadron.eta()
                  << " phi " << hadron.phi()
                  << " mothers (" << hadron.mother1()
                  << ", " << hadron.mother2() << ")"
                  << std::endl;

        const std::vector<int> contributors = hadron.motherList();
        for (int iMo : contributors) {
            if (!validEventIndex(event, iMo)) continue;

            const Particle& mo = event[iMo];
            const double dR = deltaRTruth(stable, mo);

            std::vector<int> promptAncestors;
            std::set<int> visited;
            getPromptAncestorsRecursiveImpl(event, iMo, ttbarPartons,
                                            promptAncestors, visited);
            std::sort(promptAncestors.begin(), promptAncestors.end());
            promptAncestors.erase(
                std::unique(promptAncestors.begin(), promptAncestors.end()),
                promptAncestors.end());

            std::cout << "          contributor idx " << iMo
                      << " id " << mo.id()
                      << " status " << mo.status()
                      << " eta " << mo.eta()
                      << " phi " << mo.phi()
                      << " dR(stable,parent) " << dR
                      << " TTbarJets {";

            bool first = true;
            for (int ancestor : promptAncestors) {
                auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(),
                                    ancestor);
                if (it == ttbarPartons.end()) continue;

                if (!first) std::cout << ",";
                first = false;
                std::cout << std::distance(ttbarPartons.begin(), it);
            }
            std::cout << "}" << std::endl;
        }
    }

    double selectedDR = std::numeric_limits<double>::infinity();
    const int selectedParent =
        closestHadronizationParent(event, stableIndex, primaryHadrons,
                                   &selectedDR);

    if (!validEventIndex(event, selectedParent)) {
        std::cout << "        SELECTED: no usable contributor -> rejected"
                  << std::endl;
        return;
    }

    std::vector<int> selectedPromptAncestors;
    std::set<int> selectedVisited;
    getPromptAncestorsRecursiveImpl(event, selectedParent, ttbarPartons,
                                    selectedPromptAncestors, selectedVisited);
    std::sort(selectedPromptAncestors.begin(), selectedPromptAncestors.end());
    selectedPromptAncestors.erase(
        std::unique(selectedPromptAncestors.begin(),
                    selectedPromptAncestors.end()),
        selectedPromptAncestors.end());

    bool selectedAmbiguous = false;
    const int selectedJet =
        ttbarJetFromSelectedParent(event, stableIndex, selectedParent,
                                   ttbarPartons, &selectedAmbiguous);

    const Particle& selected = event[selectedParent];
    std::cout << "        SELECTED contributor idx " << selectedParent
              << " id " << selected.id()
              << " status " << selected.status()
              << " eta " << selected.eta()
              << " phi " << selected.phi()
              << " dR " << selectedDR
              << " TTbarJets {";

    bool first = true;
    for (int ancestor : selectedPromptAncestors) {
        auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(),
                            ancestor);
        if (it == ttbarPartons.end()) continue;
        if (!first) std::cout << ",";
        first = false;
        std::cout << std::distance(ttbarPartons.begin(), it);
    }

    std::cout << "} -> ";
    if (selectedJet < 0) {
        std::cout << "REJECTED (selected contributor has no TTbar assignment)";
    } else if (selectedDR > 1.0) {
        std::cout << "REJECTED (DeltaR to selected contributor > 1.0)";
    } else {
        std::cout << "assigned TTbar jet " << selectedJet;
        if (selectedAmbiguous)
            std::cout << " after TTbar ambiguity resolution";
    }
    std::cout << std::endl;
}


void printBadTopCollectionDiagnostic(Event& event,
                                     const std::vector<int>& collection,
                                     int hardBIndex,
                                     const char* label,
                                     double reconstructedMass,
                                     double pythiaMass) {
    if (!validEventIndex(event, hardBIndex)) return;

    const std::vector<int> ttbarPartons = getPartonIndicesFromTTbar(event);

    double sumPx = 0.0;
    double sumPy = 0.0;
    double sumPz = 0.0;
    double sumE  = 0.0;

    std::vector<int> constituents;
    constituents.reserve(collection.size());

    for (int idx : collection) {
        if (!validEventIndex(event, idx)) continue;
        const auto& p = event[idx];
        if (!p.isFinal()) continue;

        sumPx += p.px();
        sumPy += p.py();
        sumPz += p.pz();
        sumE  += p.e();
        constituents.push_back(idx);
    }

    const double m2 = sumE*sumE - sumPx*sumPx - sumPy*sumPy - sumPz*sumPz;
    const double collectionMass = std::sqrt(std::max(0.0, m2));

    std::sort(constituents.begin(), constituents.end(),
              [&event](int a, int b) { return event[a].e() > event[b].e(); });

    std::cout << "\n========== BAD " << label << " MASS DIAGNOSTIC ==========" << std::endl;
    std::cout << "Pythia mass: " << pythiaMass
              << ", reconstructed mass: " << reconstructedMass << std::endl;
    std::cout << "hard b index: " << hardBIndex
              << " id: " << event[hardBIndex].id()
              << " pT: " << event[hardBIndex].pT()
              << " eta: " << event[hardBIndex].eta()
              << " phi: " << event[hardBIndex].phi() << std::endl;
    std::cout << "b collection size: " << constituents.size()
              << ", b-collection-only invariant mass from Pythia px/py/pz/E: "
              << collectionMass << std::endl;

    const std::size_t nToPrint = std::min<std::size_t>(10, constituents.size());
    std::cout << "Top " << nToPrint << " b-collection constituents by energy:" << std::endl;

    for (std::size_t rank = 0; rank < nToPrint; ++rank) {
        const int idx = constituents[rank];
        const auto& p = event[idx];
        const double dR = deltaRTruth(p, event[hardBIndex]);

        std::cout << "  [" << rank << "] particle " << idx
                  << " id " << p.id()
                  << " status " << p.status()
                  << " pT " << p.pT()
                  << " eta " << p.eta()
                  << " phi " << p.phi()
                  << " E " << p.e()
                  << " dR(b) " << dR
                  << " mothers (" << p.mother1() << ", " << p.mother2() << ")"
                  << std::endl;

        const std::vector<int> path = shortestMotherPathToTarget(event, idx, hardBIndex);
        if (path.empty()) {
            std::cout << "      NO MOTHER PATH TO HARD b FOUND" << std::endl;
        } else {
            std::cout << "      shortest mother path (stable -> hard b):" << std::endl;
            for (int pathIdx : path) {
                const auto& q = event[pathIdx];
                std::cout << "        idx " << pathIdx
                          << " id " << q.id()
                          << " status " << q.status()
                          << " mothers (" << q.mother1() << ", " << q.mother2() << ")"
                          << " pT " << q.pT()
                          << " eta " << q.eta()
                          << " phi " << q.phi()
                          << " E " << q.e()
                          << std::endl;
            }
        }

        if (ttbarPartons.size() == 6) {
            printHadronizationParentDecisionDiagnostic(
                event, idx, ttbarPartons);
        }
    }

    std::cout << "====================================================\n" << std::endl;
}


struct GenClusterTTbarClassification {
    int genClusterIndex = -1;
    std::vector<int> finalPartonAncestors;
    std::vector<int> ttbarJetIndices;      // entries 0..5
    std::vector<int> ttbarPartonIndices;   // Pythia indices for those truth quarks

    bool isMatched() const { return !ttbarJetIndices.empty(); }
    bool isMixed() const { return ttbarJetIndices.size() > 1; }
};

std::vector<GenClusterTTbarClassification> classifyGenClustersByTTbarPartonAncestry(
    Event& event,
    const vector<vector<int>>& genClusters) {

    const std::vector<int> ttbarPartons = getPartonIndicesFromTTbar(event);
    std::vector<GenClusterTTbarClassification> classifications;
    classifications.reserve(genClusters.size());

    if (ttbarPartons.size() != 6)
        return classifications;

    for (std::size_t g = 0; g < genClusters.size(); ++g) {
        GenClusterTTbarClassification cls;
        cls.genClusterIndex = static_cast<int>(g);

        if (!genClusters[g].empty()) {
            getAncestors(event, genClusters[g][0], cls.finalPartonAncestors);
        }

        std::set<int> matchedTruthJets;

        for (int finalParton : cls.finalPartonAncestors) {
            std::vector<int> promptAncestors;
            std::set<int> visited;
            getPromptAncestorsRecursiveImpl(event, finalParton, ttbarPartons,
                                            promptAncestors, visited);

            std::sort(promptAncestors.begin(), promptAncestors.end());
            promptAncestors.erase(
                std::unique(promptAncestors.begin(), promptAncestors.end()),
                promptAncestors.end());

            for (int ancestor : promptAncestors) {
                auto it = std::find(ttbarPartons.begin(), ttbarPartons.end(), ancestor);
                if (it == ttbarPartons.end()) continue;
                matchedTruthJets.insert(
                    static_cast<int>(std::distance(ttbarPartons.begin(), it)));
            }
        }

        for (int truthJet : matchedTruthJets) {
            cls.ttbarJetIndices.push_back(truthJet);
            cls.ttbarPartonIndices.push_back(ttbarPartons[truthJet]);
        }

        classifications.push_back(cls);
    }

    return classifications;
}

void printGenClusterTTbarClassifications(
    Event& event,
    const vector<vector<int>>& genClusters) {

    const auto classifications =
        classifyGenClustersByTTbarPartonAncestry(event, genClusters);

    std::cout << "genCluster -> TTbar truth-jet classification:" << std::endl;

    for (const auto& cls : classifications) {
        std::cout << "  genCluster " << cls.genClusterIndex << ": final-parton ancestors {";
        for (std::size_t i = 0; i < cls.finalPartonAncestors.size(); ++i) {
            const int idx = cls.finalPartonAncestors[i];
            if (i) std::cout << ", ";
            std::cout << idx << "(id " << event[idx].id() << ")";
        }
        std::cout << "} -> ";

        if (!cls.isMatched()) {
            std::cout << "no TTbar jet";
        } else {
            std::cout << "TTbar jet" << (cls.ttbarJetIndices.size() == 1 ? " " : "s {");
            for (std::size_t i = 0; i < cls.ttbarJetIndices.size(); ++i) {
                if (i) std::cout << ", ";
                const int jet = cls.ttbarJetIndices[i];
                const int parton = cls.ttbarPartonIndices[i];
                std::cout << jet << " [parton " << parton
                          << ", id " << event[parton].id() << "]";
            }
            if (cls.ttbarJetIndices.size() > 1) std::cout << "} [MIXED]";
        }

        std::cout << std::endl;
    }

    // Also print the inverse map: where each of the six TTbar jets appears.
    const std::vector<int> ttbarPartons = getPartonIndicesFromTTbar(event);
    if (ttbarPartons.size() != 6) return;

    std::cout << "TTbar truth jet -> genCluster classification:" << std::endl;
    for (int t = 0; t < 6; ++t) {
        std::vector<int> matchedClusters;
        for (const auto& cls : classifications) {
            if (std::find(cls.ttbarJetIndices.begin(), cls.ttbarJetIndices.end(), t)
                != cls.ttbarJetIndices.end()) {
                matchedClusters.push_back(cls.genClusterIndex);
            }
        }

        std::cout << "  TTbar jet " << t
                  << " parton " << ttbarPartons[t]
                  << " (id " << event[ttbarPartons[t]].id() << ") -> ";

        if (matchedClusters.empty()) {
            std::cout << "no genCluster";
        } else {
            for (std::size_t i = 0; i < matchedClusters.size(); ++i) {
                if (i) std::cout << ", ";
                std::cout << "genCluster " << matchedClusters[i];
            }
        }
        std::cout << std::endl;
    }
}

struct TTbarClusterMatch {
    int ttbarJetIndex = -1;       // 0..5 in getPartonIndicesFromTTbar ordering
    int ttbarPartonIndex = -1;    // Pythia event-record index of truth quark
    int genClusterIndex = -1;     // best-overlap genCluster
    unsigned sharedParticles = 0;
    double overlapFraction = 0.0; // shared / particles in TTbar truth cluster
};

std::vector<TTbarClusterMatch> matchTTbarTruthClustersToGenClusters(
    Event& event,
    const vector<vector<int>>& genClusters,
    const vector<vector<int>>& ttbarClusters) {

    const std::vector<int> ttbarPartons = getPartonIndicesFromTTbar(event);
    std::vector<TTbarClusterMatch> matches;

    if (ttbarPartons.size() != 6 || ttbarClusters.size() != 6)
        return matches;

    matches.reserve(6);

    for (std::size_t t = 0; t < 6; ++t) {
        TTbarClusterMatch match;
        match.ttbarJetIndex = static_cast<int>(t);
        match.ttbarPartonIndex = ttbarPartons[t];

        std::set<int> truthParticles(ttbarClusters[t].begin(),
                                     ttbarClusters[t].end());

        for (std::size_t g = 0; g < genClusters.size(); ++g) {
            unsigned shared = 0;
            for (int particle : genClusters[g]) {
                if (truthParticles.count(particle)) ++shared;
            }

            if (shared > match.sharedParticles) {
                match.sharedParticles = shared;
                match.genClusterIndex = static_cast<int>(g);
            }
        }

        if (!truthParticles.empty()) {
            match.overlapFraction =
                static_cast<double>(match.sharedParticles) /
                static_cast<double>(truthParticles.size());
        }

        matches.push_back(match);
    }

    return matches;
}

void printTTbarGenClusterMatches(Event& event,
                                 const vector<vector<int>>& genClusters,
                                 const vector<vector<int>>& ttbarClusters) {
    const auto matches = matchTTbarTruthClustersToGenClusters(
        event, genClusters, ttbarClusters);

    std::cout << "TTbar -> genCluster overlap matches:" << std::endl;
    for (const auto& match : matches) {
        std::cout << "  TTbar jet " << match.ttbarJetIndex
                  << " parton " << match.ttbarPartonIndex
                  << " (id " << event[match.ttbarPartonIndex].id() << ")"
                  << " -> genCluster " << match.genClusterIndex
                  << ", shared particles " << match.sharedParticles
                  << ", truth overlap " << match.overlapFraction
                  << std::endl;
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


void getStableDescendantsRecursive(Event& event, int iParton, std::vector<int>& products) {
    if (iParton <= 0 || iParton >= event.size()) return;

    const auto& p = event[iParton];
    if (p.isFinal()) {
        products.push_back(iParton);
        return;
    }

    const int firstDaughter = p.daughter1();
    const int lastDaughter  = p.daughter2();
    if (firstDaughter <= 0 || lastDaughter <= 0) return;

    const int first = std::max(1, std::min(firstDaughter, lastDaughter));
    const int last  = std::min(event.size() - 1, std::max(firstDaughter, lastDaughter));

    for (int iDau = first; iDau <= last; ++iDau)
        getStableDescendantsRecursive(event, iDau, products);
}
void getStableDescendants(Event& event, int iParton,std::vector<int>& products) {
    getStableDescendantsRecursive(event, iParton, products);
    std::sort(products.begin(),products.end());
    auto it_p = std::unique(products.begin(),products.end());
    products.erase(it_p, products.end());
}

void getFinalPartonDescendantsRecursive(Event& event, int iParton, std::vector<int>& products) {
    if (iParton <= 0 || iParton >= event.size()) return;

    const auto& p = event[iParton];
    if (p.isFinalPartonLevel()) {
        products.push_back(iParton);
        return;
    }

    const int firstDaughter = p.daughter1();
    const int lastDaughter  = p.daughter2();
    if (firstDaughter <= 0 || lastDaughter <= 0) return;

    const int first = std::max(1, std::min(firstDaughter, lastDaughter));
    const int last  = std::min(event.size() - 1, std::max(firstDaughter, lastDaughter));

    for (int iDau = first; iDau <= last; ++iDau)
        getFinalPartonDescendantsRecursive(event, iDau, products);
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

	evt.invMasses.clear();
	evt.recMasses.clear();

	// Get the decaying top and anti-top.
	const int iTop     = findLastParticleWithId(*evt.pythiaEvent,  6);
	const int iAntiTop = findLastParticleWithId(*evt.pythiaEvent, -6);

	int b = -1;
	int wPlus = -1;
	int bbar = -1;
	int wMinus = -1;

	if (iTop >= 0 && iAntiTop >= 0 &&
    		findTopDecay(*evt.pythiaEvent, iTop, b, wPlus) &&
    		findTopDecay(*evt.pythiaEvent, iAntiTop, bbar, wMinus)) {

    		std::vector<double> masses;

    		evt.invMasses.push_back((*evt.pythiaEvent)[wPlus].m());
    		evt.invMasses.push_back((*evt.pythiaEvent)[wMinus].m());
    		evt.invMasses.push_back((*evt.pythiaEvent)[iTop].m());
    		evt.invMasses.push_back((*evt.pythiaEvent)[iAntiTop].m());
	}


	std::vector<std::size_t> rejectedMixedPerTTbarJet;
	getFinalParticleClustersFromTTbar(*evt.pythiaEvent,
	                                  evt.genClustersFromHardCollision,
	                                  evt.hardScatterSize,
	                                  &rejectedMixedPerTTbarJet);

	// Keep the legacy evt.genClusters interface, but make it represent the
	// six fixed TTbar truth jets.  From this point onward:
	//
	//   genClusters[0] = b
	//   genClusters[1] = W daughter
	//   genClusters[2] = W daughter
	//   genClusters[3] = bbar
	//   genClusters[4] = W daughter
	//   genClusters[5] = W daughter
	//
	// This preserves older code that reads evt.genClusters while making its
	// indexing identical to evt.genClustersFromHardCollision and evt.genJets.
	evt.genClusters = evt.genClustersFromHardCollision;

	// Reconstruct W+/W-/top/anti-top masses from the six TTbar truth
	// particle collections.  Keep the same ordering as evt.invMasses:
	//   0 = W+, 1 = W-, 2 = top, 3 = anti-top.
	if (evt.genClustersFromHardCollision.size() >= 6) {
		std::vector<rk::P4> ttbarJetP4s(6);

		for (std::size_t i = 0; i < 6; ++i) {
			for (int pid : evt.genClustersFromHardCollision[i]) {
				if (pid <= 0 || pid >= evt.pythiaEvent->size()) continue;

				const Pythia8::Particle& p = (*evt.pythiaEvent)[pid];
				if (!p.isFinal()) continue;

				rk::P4 part = rk::P4(
					p.pT() * geom3::Vector3(
						std::cos(p.phi()),
						std::sin(p.phi()),
						std::sinh(p.eta())),
					p.m());

				ttbarJetP4s[i] = ttbarJetP4s[i] + part;
			}
		}

		const rk::P4 wPlusRec = ttbarJetP4s[1] + ttbarJetP4s[2];
		const rk::P4 wMinusRec = ttbarJetP4s[4] + ttbarJetP4s[5];
		const rk::P4 topRec = ttbarJetP4s[0] + wPlusRec;
		const rk::P4 antiTopRec = ttbarJetP4s[3] + wMinusRec;

		evt.recMasses.push_back(wPlusRec.m());
		evt.recMasses.push_back(wMinusRec.m());
		evt.recMasses.push_back(topRec.m());
		evt.recMasses.push_back(antiTopRec.m());

		// Compact one-line diagnostics for the two b collections.  These are
		// deliberately printed for every event so they can be parsed into 2D
		// correlations with reconstructed top mass.
		const std::size_t rejectedB =
		    rejectedMixedPerTTbarJet.size() > 0 ? rejectedMixedPerTTbarJet[0] : 0;
		const std::size_t rejectedBbar =
		    rejectedMixedPerTTbarJet.size() > 3 ? rejectedMixedPerTTbarJet[3] : 0;
		printBCollectionSummary(*evt.pythiaEvent,
		                        evt.genClustersFromHardCollision[0],
		                        b, "B", topRec.m(), rejectedB);
		printBCollectionSummary(*evt.pythiaEvent,
		                        evt.genClustersFromHardCollision[3],
		                        bbar, "BBAR", antiTopRec.m(), rejectedBbar);

		// Diagnose the long reconstructed top-mass tails without changing any particle assignment.
		constexpr double badTopMassThreshold = 250.0;
		if (topRec.m() > badTopMassThreshold && evt.invMasses.size() >= 4) {
			printBadTopCollectionDiagnostic(*evt.pythiaEvent,
			                                evt.genClustersFromHardCollision[0],
			                                b,
			                                "TOP",
			                                topRec.m(),
			                                evt.invMasses[2]);
		}
		if (antiTopRec.m() > badTopMassThreshold && evt.invMasses.size() >= 4) {
			printBadTopCollectionDiagnostic(*evt.pythiaEvent,
			                                evt.genClustersFromHardCollision[3],
			                                bbar,
			                                "ANTITOP",
			                                antiTopRec.m(),
			                                evt.invMasses[3]);
		}
	}
	printGenClusterTTbarClassifications(*evt.pythiaEvent, evt.genClusters);
	printTTbarGenClusterMatches(*evt.pythiaEvent, evt.genClusters, evt.genClustersFromHardCollision);
        evt.genJetsReady = true;

        std::vector<rk::P4> finalParticles;
        finalParticles.reserve(evt.pythiaEvent->size());

        // Fill genEvent from the six TTbar truth-particle collections rather
        // than from the generic genClusters.  This preserves the same event
        // representation used previously, but the input particles are now
        // exactly the particles assigned to:
        //   0 = b, 1 = W q, 2 = W q,
        //   3 = bbar, 4 = W q, 5 = W q.
        //
        // Build evt.genJets from the same six collections so evt.genJets[0..5]
        // has the identical fixed TTbar ordering.
	for (long unsigned int i=0; i<evt.genClusters.size(); ++i) {
		rk::P4 genJet;
                if (!evt.genClusters[i].empty()) {
			int jetParts = evt.genClusters[i].size();
                        for (int j=0; j<jetParts; ++j) {
				int pid = evt.genClusters[i][j];
				const Pythia8::Particle& p = (*evt.pythiaEvent)[pid];
				if (p.isFinal()) {
					rk::P4 part = rk::P4(
                                        p.pT()*geom3::Vector3(
                                            cos(p.phi()),
                                            sin(p.phi()),
                                            sinh(p.eta())),
                                        p.m());
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
