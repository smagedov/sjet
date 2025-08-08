set outdir = JetHistories
mkdir -p $outdir
# set dumpFile = ${outdir}/jetEnergyFlow.txt
set outputPrefix = ${outdir}/jethist
set avePtRef = 10.0
set muRef = 1000
set avePtPileup = 2.0
set muPileup = 0.0
set nEvents = 100
set distanceCutoff = 100.0
set seed = 115371470657
set countPeriod = 10

./gaussianJetStudy --countPeriod $countPeriod --seed $seed $avePtRef $muRef $avePtPileup $muPileup $nEvents $distanceCutoff $outputPrefix
