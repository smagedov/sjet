set avePtRef = 10.0
set muRef = 20
set avePtProbe = 10.0
set muProbe = 20
set avePtPileup = 2.0
set muPileup = 0.0
set probeJetWidth = 0.4
set deltaR = 1.5
set nEvents = 1000
set distanceCutoff = 0.4
set seed = 115371470657
set countPeriod = 10
set dumpFile = "dijetEnergyFlow.txt"
set outputFile = "dijetStudyExampleOutput.txt"
set histFile = "dijetStudyClusteringHistory.txt"

./gaussianDijetClustComp -r --countPeriod $countPeriod --seed $seed --dumpFile $dumpFile $avePtRef $muRef $avePtProbe $muProbe $avePtPileup $muPileup $probeJetWidth $deltaR $nEvents $distanceCutoff $outputFile $histFile
