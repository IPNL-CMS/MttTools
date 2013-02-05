#! /bin/sh

echo "Usage: generatePUProfileForData.sh [crab_json_file] [input_lumi_json] [output_file]"
echo
echo "Generating PU profile..."

pileupCalc.py -i $1 --inputLumiJSON $2 --calcMode=true --maxPileupBin=70 --numPileupBins=70 --minBiasXsec 69300 $3_nominal.root
pileupCalc.py -i $1 --inputLumiJSON $2 --calcMode=true --maxPileupBin=70 --numPileupBins=70 --minBiasXsec 72765 $3_puUp.root
pileupCalc.py -i $1 --inputLumiJSON $2 --calcMode=true --maxPileupBin=70 --numPileupBins=70 --minBiasXsec 65835 $3_puDown.root
