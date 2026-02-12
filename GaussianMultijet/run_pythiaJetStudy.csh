set nEvents = 1
set nJets = 12
set distanceCutoff = 0.6
set maxDistance = 10.5
set seed = 115371470657
set countPeriod = 10
set dumpFile = "pythiaEnergyFlow.txt"
set outputFile = "pythiaStudyExampleOutput.txt"
set inputFile = "TTbar.cmnd"

./pythiaJetStudy --countPeriod $countPeriod --seed $seed --dumpFile $dumpFile $nEvents $nJets $distanceCutoff $outputFile $inputFile
