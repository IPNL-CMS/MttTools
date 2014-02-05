#! /bin/bash

PLOTIT="../plotIt/plotIt"
date=`LANG=C date +%d%b%y`

# Create folder structure
for type in semimu semie; do
  for b in 1 2; do
    mkdir -p "plots/${date}/${b}-btag/${type}"
  done
done

for type in semimu semie; do
  for b in 1 2; do
    config="MTT_${b}btag_${type}.yml"
    output="plots/${date}/${b}-btag/${type}"

    $PLOTIT -o $output $config
  done
done
