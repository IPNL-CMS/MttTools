#! /bin/bash

for mass in 500 750 1000 1250 1500 2000; do
  cat MTT_Zprime_${mass}_Narrow_*semimu.list > MTT_Zprime_${mass}_Narrow_full_stat_semimu.list
  cat MTT_Zprime_${mass}_Narrow_*semie.list > MTT_Zprime_${mass}_Narrow_full_stat_semie.list
  cat MTT_Zprime_${mass}_Large_*semimu.list > MTT_Zprime_${mass}_Large_full_stat_semimu.list
  cat MTT_Zprime_${mass}_Large_*semie.list > MTT_Zprime_${mass}_Large_full_stat_semie.list
done;
