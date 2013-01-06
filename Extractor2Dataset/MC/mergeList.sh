#! /bin/bash

for mass in 750 1000 1250 1500 2000; do
  cat MTT_Zprime_${mass}_Narrow_*_semimu.list > MTT_Zprime_${mass}_Narrow_full_stat_semimu.list
  cat MTT_Zprime_${mass}_Narrow_*_semie.list > MTT_Zprime_${mass}_Narrow_full_stat_semie.list
done;
