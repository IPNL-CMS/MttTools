#! /bin/sh

echo "Usage: generatePUProfileForData.sh [input_lumi_json]"
echo
echo "Generating PU profiles..."

date=`date +%d%b%y`

# Single Mu

for a in SingleMu_Run2012*.json; do
  run=${a#Single*_};
  run=${run%.*}
  output="pileup_SingleMu_${run}_${date}_"

  for syst in nominal puUp puDown; do
    case $syst in
      nominal)
        xsec=69400
        ;;
      puUp)
        xsec=72870
        ;;
      puDown)
        xsec=65930
        ;;
    esac
    echo "Generating PU distribution for ${run}:${syst}"
    pileupCalc.py -i $a --inputLumiJSON $1 --calcMode=true --maxPileupBin=100 --numPileupBins=100 --minBiasXsec $xsec ${output}${syst}.root
    var="outputList_${syst}"
    declare "${var}=${!var} ${output}${syst}.root"
  done
done

hadd pileup_SingleMu_full_stat_${date}_nominal.root ${outputList_nominal}
hadd pileup_SingleMu_full_stat_${date}_puUp.root ${outputList_puUp}
hadd pileup_SingleMu_full_stat_${date}_puDown.root ${outputList_puDown}

for syst in nominal puUp puDown; do
  var="outputList_${syst}"
  declare "${var}=\"\""
done

for a in SingleElectron_Run2012*.json; do
  run=${a#Single*_};
  run=${run%.*}
  output="pileup_SingleElectron_${run}_${date}_"

  for syst in nominal puUp puDown; do
    case $syst in
      nominal)
        xsec=69400
        ;;
      puUp)
        xsec=72870
        ;;
      puDown)
        xsec=65930
        ;;
    esac
    echo "Generating PU distribution for ${run}:${syst}"
    pileupCalc.py -i $a --inputLumiJSON $1 --calcMode=true --maxPileupBin=100 --numPileupBins=100 --minBiasXsec $xsec ${output}${syst}.root
    var="outputList_${syst}"
    declare "${var}=${!var} ${output}${syst}.root"
  done
done

hadd pileup_SingleElectron_full_stat_${date}_nominal.root ${outputList_nominal}
hadd pileup_SingleElectron_full_stat_${date}_puUp.root ${outputList_puUp}
hadd pileup_SingleElectron_full_stat_${date}_puDown.root ${outputList_puDown}
