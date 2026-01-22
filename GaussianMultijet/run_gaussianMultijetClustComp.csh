set avePtRef = 10.0
set muRef = 20
set avePtProbe = 10.0
set muProbe = 20
set avePtPileup = 2.0
set muPileup = 0.0
set probeJetWidth = 0.4
set deltaR = 3.0
set nEvents = 10
set nJets = 3
set distanceCutoff = 0.5
set seed = 115371470657
set countPeriod = 10
set dumpFile = "dijetEnergyFlow.txt"
set outputFile = "dijetStudyExampleOutput.txt"

./gaussianMultijetClustComp -r --countPeriod $countPeriod --seed $seed --dumpFile $dumpFile $avePtRef $muRef $avePtProbe $muProbe $avePtPileup $muPileup $probeJetWidth $deltaR $nEvents $nJets $distanceCutoff $outputFile
