#! /bin/bash

PLOTIT="../plotIt/plotIt"
date=`LANG=C date +%d%b%y`

# Create folder structure
for type in semimu semie; do
  for b in 1 2; do
    mkdir -p "plots/${date}/QCD/${b}-btag/${type}"
  done
done

for type in semimu semie; do
  for b in 1 2; do
    config="MTT_${b}btag_${type}_QCD.yml"
    output="plots/${date}/QCD/${b}-btag/${type}"

    if [ -f $config ]; then
      $PLOTIT -o $output $config
    fi
  done
done
