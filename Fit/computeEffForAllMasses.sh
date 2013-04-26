#! /bin/bash

for syst in "" "--jec up" "--jec down" "--jer up" "--jer down" "--pu up" "--pu down" "--pdf up" "--pdf down"; do
  ./computeEff --b-tag 2 $syst
  ./computeEff --b-tag 1 $syst
done;
