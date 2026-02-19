This directory holds the necessary files for working with pythia files, and
running stablejet clustering on the finalized particles outputted by pythia.
The main goal of these files is to covert the information obtained from
pythia into a form useable by the stablejet clustering.

In order to run the analysis framework, one needs to run the corresponding .csh
file.

csh ./run_pythiaJetStudy.csh

Within run_pythiaJetStudy.csh there are several options which can be used to control
the clustering.

The main body of the analysis code is found in pythiaJetStudy.cc which calls several
modules to perform analysis.

PythiaEventMaker	Utilizes the given pythia card to generate a pythia event
			and appends it under evt.PythiaEvent as a member of the
			"Event" to be analyzed.

PythiaEventTTbarGen	Utilizing the pythia event generated in PythiaEventMaker,
			populates the relevant members of the "Event" utilized in
			clustering. Utilizing information from the derived
			genClusters, it generates the necessary genEvent and genJets
			information.

RunDiffusionClustering	Runs the "slow" clustering algorithm on the generated
			genEvent.

PythiaJetMaker		Having run the clustering, it then creates the appropriate
			clustered jets and stores them in simpleDiffusionJets member
			of the "Event".

PythiaOutputPrinter	Prints out the information about the final state particles
			of the pythia event, along with which genJets the particles
			are assigned to.

GenParticleDump		This module can be used to dump the
			particles created by DijetMaker into a file
			(in a text form). This is useful for
			visualizing the average energy flow.

FullJetHistoryPrinter	Outputs the full clustering history for a particular clustred 
			jet, and writes it to a designated .txt file.

PythiaClustering	Generates all possible permutations of jet matching,
			and minimizes the dR for the overall event. Storing only the
			final, minimized jet matching.

JetHistoryCopy		Creates a copy of the clustering history, replacing each
			rk::P4 object with a vector that seperates the contribution
			of each original genJet for each cluster in the history.
