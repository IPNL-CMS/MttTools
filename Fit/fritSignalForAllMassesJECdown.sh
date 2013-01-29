#! /bin/bash

for mass in 500 750 1000 1250 1500; do
  echo "Fit signal for M=$mass"
  ./fritSignal -m $mass --b-tag 2 -i ../Extractor2Dataset/MTT_Signal_Zprime_${mass}_Narrow_*_merged_*_JECdown.root --jec JECdown
  ./fritSignal -m $mass --b-tag 1 -i ../Extractor2Dataset/MTT_Signal_Zprime_${mass}_Narrow_*_merged_*_JECdown.root --jec JECdown
done;
