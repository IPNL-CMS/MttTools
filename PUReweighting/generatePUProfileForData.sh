#! /bin/sh

echo "Usage: generatePUProfileForData.sh [crab_json_file] [input_lumi_json] [output_file]"
echo
echo "Generating PU profile..."

pileupCalc.py -i $1 --inputLumiJSON $2 --calcMode=true --maxPileupBin=70 --numPileupBins=70 --minBiasXsec 69300 $3
#pileupCalc.py -i $1 --inputLumiJSON $2 --calcMode=true --maxPileupBin=70 --numPileupBins=70 --minBiasXsec 70000 $3
