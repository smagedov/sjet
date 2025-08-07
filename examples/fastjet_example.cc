#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <algorithm>

#include "fastjet/config.h"
#include "sjet/DiffusionDistFastJet.hh"

using namespace std;

int main(){
//	unsigned nParticles;
//	string filename = "single-event.dat";
//        ifstream inputFile(filename);
//
//        vector<vector<double>> particleArray;
//
//        string line;
//        while (getline(inputFile, line)) {
//                vector<double> particle;
//                stringstream ss(line);
//                double value;
//                while (ss >> value) {
//                        particle.push_back(value);
//                }
//                particleArray.push_back(particle);
//        }
//
//        inputFile.close();
//
//       nParticles = particleArray.size();
//       cout << "Total number of particles found: " << nParticles << endl;
//
	vector<fastjet::PseudoJet> input_particles;
//
//	for (unsigned i = 0; i < nParticles; ++i) {
//		input_particles.push_back(fastjet::PseudoJet(particleArray[i][0], particleArray[i][1], particleArray[i][2], particleArray[i][3]));
//	}

	double px, py , pz, E;
	while (cin >> px >> py >> pz >> E) {
		input_particles.push_back(fastjet::PseudoJet(px,py,pz,E));
	}

	cout << input_particles.size() << endl;

	double maxDist = 0.25;
	sjet::DiffusionDistFastJet sjet(maxDist);
	fastjet::JetDefinition jet_def(&sjet);

	fastjet::ClusterSequence clust_seq(input_particles, jet_def);

	double ptmin = 0.0;
	vector<fastjet::PseudoJet> exclusive_jets = clust_seq.exclusive_jets(maxDist);

	cout << "Ran " << jet_def.description() << endl;

	cout << "\n\nTotal number of jets obtained: " << clust_seq.n_exclusive_jets(maxDist) << endl;

	printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
	 
	for (unsigned int i = 0; i < exclusive_jets.size(); i++) {
		printf("%5u %15.8f %15.8f %15.8f\n",
				i, exclusive_jets[i].rap(), exclusive_jets[i].phi(),
				exclusive_jets[i].perp());
	}

	return 0;

}
