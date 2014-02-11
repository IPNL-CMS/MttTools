#! /bin/bash

for mass in 400 500 600 700 800; do
  echo "Fit signal for M=$mass"
  ./higgsSignalToPDF -m $mass --b-tag 2 -i ../Extractor2Dataset/datasets/Latest/Signal_S0_S_i_M${mass}_cpl1_pseudoscalar_histos_nominal_merged.root
  ./higgsSignalToPDF -m $mass --b-tag 1 -i ../Extractor2Dataset/datasets/Latest/Signal_S0_S_i_M${mass}_cpl1_pseudoscalar_histos_nominal_merged.root
done;
