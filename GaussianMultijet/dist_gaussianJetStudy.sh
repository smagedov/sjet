#!/bin/bash

for f in $(seq 0.1 0.1 3.0); do

	mkdir dijetDisStudy/$f

	sed -i "s/^set distanceCutoff = .*/set distanceCutoff = $f/" run_gaussianDijetStudy.csh

	for i in $(seq 0.1 0.1 3.0); do
		echo "$i"

		sed -i "s/^set deltaR = .*/set deltaR = $i/" run_gaussianDijetStudy.csh

		csh run_gaussianDijetStudy.csh

		mv dijetEnergyFlow.txt "dijetEnergyFlow_${i}.txt"
		mv dijetStudyExampleOutput.txt "dijetStudyExampleOutput_${i}.txt"

		mv *.txt dijetDisStudy/$f/
	done
done
