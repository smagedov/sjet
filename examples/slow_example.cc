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

#include "CmdLine.hh"

using namespace std;
using namespace sjet;

typedef ClusteringSequence<DistanceCalculator, rk::P4> MyClusteringSequence;

static void print_usage(const char* progname) {
        cout << "\nUsage: " << progname << " maxDistance\n" << endl;
}

rk::P4 randomParticle(const double px, const double py, const double pz, const double E) {

//	const double pt = drand48()*9 + 1;
//	const double eta = drand48()*4 - 2;
//	const double phi = drand48()*M_PI*2;
//	const double mass = drand48()*4 + 1;
//	return rk::P4(pt*geom3::Vector3(cos(phi), sin(phi), sinh(eta)), mass);

	return rk::P4(geom3::Vector3(px, py, pz), E);
}

static void testing(const unsigned nParticles, const double maxDistance, vector<vector<double>>& particleArray) {
	vector<rk::P4> particles;
	particles.reserve(nParticles);
	for (unsigned i = 0; i < nParticles; ++i) {
		particles.push_back(randomParticle(particleArray[i][0], particleArray[i][1], particleArray[i][2], particleArray[i][3]));
	}
	DistanceCalculator calc;
	MyClusteringSequence clustSequence(calc);
	clustSequence.init(particles);
	if (maxDistance > 0.0) {
		clustSequence.run(maxDistance);
	} else {
		clustSequence.run();
	}

	cout << "Max d: " << maxDistance << " Final nClusters: " << clustSequence.nClusters() << endl;
}

int main(int argc, char *argv[]) {
	CmdLine cmdline(argc, argv);
	
	if (argc == 1) {
		print_usage(cmdline.progname());
		return 0;
	}

	unsigned nParticles;
	double maxDistance;

	try {
		cmdline.optend();

		if (cmdline.argc() != 1)
			throw CmdLineError("wrong number of command line arguments");
		cmdline >> maxDistance;

		if (maxDistance < 0.0)
			throw CmdLineError("Maximum distance must be non-negative");
	}

	catch (const CmdLineError& e) {
		std::cerr << "Error in " << cmdline.progname() << ": "
			<< e.str() << std::endl;
		return 1;
	}

        string filename = "single-event.dat";
        ifstream inputFile(filename);

        if (!inputFile.is_open()) {
                throw CmdLineError("Error: Could not open file ");
        }

        vector<vector<double>> particleArray;

        string line;
        while (getline(inputFile, line)) {
                vector<double> particle;
                stringstream ss(line);
                double value;
                while (ss >> value) {
                        particle.push_back(value);
                }
                particleArray.push_back(particle);
        }

        inputFile.close();

        nParticles = particleArray.size();
        cout << nParticles << endl;

    
 	testing(nParticles, maxDistance, particleArray);
	
	return 0;
}
