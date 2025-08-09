#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cmath>

#include "rk/rk.hh"

#include "sjet/DistanceCalculator.hh"
#include "sjet/ClusteringSequence.hh"

#include "Pythia8/Pythia.h"

#include "CmdLine.hh"

//======================================================================
// Hard-coded settings

// Jet and hadron pT thresholds.
// Will only show particles with pT > pTmin and |y| < yMax.
double pTmin_jet = 25;
double pTmin_hadron = 1;
double yMax = 4;

// Amount of pileup. Average number of inelastic pp collisions per event
// (=bunch-crossing). Set to zero to turn off pileup.
double mu = 60;

using namespace Pythia8;
using namespace std;
using namespace sjet;
// =====================================================================

typedef ClusteringSequence<DistanceCalculator, rk::P4> MyClusteringSequence;

static void print_usage(const char* progname) {
        cout << "\nUsage: " << progname << " maxDistance nEvents\n" << endl;
}

static void testing(const double maxDistance, vector<rk::P4>& particles) {
        DistanceCalculator calc;
        MyClusteringSequence clustSequence(calc);
        clustSequence.init(particles);
        if (maxDistance > 0.0) {
                clustSequence.run(maxDistance);
        } else {
                clustSequence.run();
        }

	// int hist = clustSequence.clustHist().size()-1;
	// int nClust = clustSequence.nClusters();

        // cout << "Max d: " << maxDistance << " Final nClusters: " << nClust << " History Size: " << hist << endl;
	// double stab;
	// for (int i = 0; i < nClust; i++) {
	//	stab = clustSequence.stability(hist - i);
		// cout << "stability for branch " << i << ": " << stab << endl;
	//}
}

int main(int argc, char *argv[]) {	
	CmdLine cmdline(argc, argv);

        if (argc == 1) {
                print_usage(cmdline.progname());
                return 0;
        }

        int nEvents;
        double maxDistance;

        try {
                cmdline.optend();

                if (cmdline.argc() != 2)
                        throw CmdLineError("wrong number of command line arguments");
                cmdline >> maxDistance;
		cmdline >> nEvents;

                if (maxDistance < 0.0)
                        throw CmdLineError("Maximum distance must be non-negative");
        }

        catch (const CmdLineError& e) {
                std::cerr << "Error in " << cmdline.progname() << ": "
                        << e.str() << std::endl;
                return 1;
        }

	Pythia pythia;
	pythia.readString("Beams:eCM = 13600.");
	pythia.readString("HiggsSM:ffbar2HW = on");
	// Force H->bb decays and hadronic W decays.
	pythia.readString("25:onMode = off");
	pythia.readString("25:onIfAny = 5");
	pythia.readString("24:onMode = off");
	pythia.readString("24:onIfAny = 1 2 3 4 5");

	// If Pythia fails to initialize, exit with error.
	if (!pythia.init()) return 1;

	// Pileup particles
	Pythia pythiaPU;
	pythiaPU.readString("Beams:eCM = 13600.");
	pythiaPU.readString("SoftQCD:inelastic = on");
	if (mu > 0) pythiaPU.init();

	std::vector<std::vector<rk::P4>> stbl_events;

	auto &event = pythia.event;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
		if (!pythia.next()) continue;
		
		// Identify particles. Jets are built from all stable particles after
		// hadronization (particle-level jets).
		std::vector<Particle> VH, ptcls_hs, ptcls_pu;
		std::vector<rk::P4> stbl_ptcls;
		for (int i = 0; i < event.size(); ++i) {
			auto &p = event[i];
			if (p.isResonance() && p.status() == -62) VH.push_back(p);
			if (not p.isFinal()) continue;
			stbl_ptcls.push_back(rk::P4(geom3::Vector3(p.px(), p.py(), p.pz()), p.e()));
			ptcls_hs.push_back(p);
		}

		stbl_events.push_back(stbl_ptcls);
		
		// Should not happen!
		if (VH.size()!= 2) continue;
	}
		
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;

	auto t1 = high_resolution_clock::now();
	for (int i = 0; i < stbl_events.size(); ++i) {
		testing(maxDistance, stbl_events[i]);
	}
	auto t2 = high_resolution_clock::now();

	duration<double, std::milli> ms_double = t2 - t1;
	std::cout << ms_double.count()/nEvents << " ms" << std::endl;

	return 0;
}
