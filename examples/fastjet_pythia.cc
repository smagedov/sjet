#include <chrono>
#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/config.h"
#include "sjet/DiffusionDistFastJet.hh"
#include "CmdLine.hh"

//==========================================================================

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

//==========================================================================

static void print_usage(const char* progname) {
        cout << "\nUsage: " << progname << " maxDistance nEvents\n" << endl;
}

// Example main program to vizualize jet algorithms.

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
	// Description of the process (using ROOT's TLatex notation).
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

	stab::DiffusionDistFastJet stab(maxDistance);
        fastjet::JetDefinition jet_def(&stab);;

	std::vector<std::vector<fastjet::PseudoJet>> stbl_events;

	auto &event = pythia.event;
	for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
		if (!pythia.next()) continue;
		
		// Identify particles. Jets are built from all stable particles after
		// hadronization (particle-level jets).
		std::vector<Particle> VH, ptcls_hs, ptcls_pul;
		std::vector<fastjet::PseudoJet> stbl_ptcls;
		for (int i = 0; i < event.size(); ++i) {
			auto &p = event[i];
			if (p.isResonance() && p.status() == -62) VH.push_back(p);
			if (not p.isFinal()) continue;
			stbl_ptcls.push_back(fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.e()));
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
		fastjet::ClusterSequence clustSeq(stbl_events[i], jet_def);
	}
	auto t2 = high_resolution_clock::now();

	duration<double, std::milli> ms_double = t2 - t1;
	std::cout << ms_double.count()/nEvents << " ms" << std::endl;
	
	return 0;
}
