#! /bin/bash

echo "Are you sure you want to clean the current directory? All ROOT files, images and text files will be removed. Continue? [y/n]"
read ok

if [ "$ok" == "y" ]; then
  find . -maxdepth 1 -type f -name '*.json' -not -name 'parameters.json' | xargs rm
  rm *.lock;
  rm *.pdf;
  rm *.png;
  rm *.gif;
  rm *.jpg;
  find . -maxdepth 1 -type f -name '*.txt' -not -name 'ds_*.txt' | xargs rm
fi
