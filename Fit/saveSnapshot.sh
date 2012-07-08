#! /bin/bash

# This script save a snapshot of the current analysis in a folder named with the current `date`

if [ -z "$1" ]; then
  folder="snapshots/`date +%Y-%m-%d@%H:%M:%S`"
else
  folder="snapshots/$1"
fi

if [ -d "$folder" ]; then
  echo "Error: a folder named $folder already exists. Exiting."
  exit 1
fi

mkdir "$folder"
mkdir -p "$folder/toys/results/"

# Start saving
echo "Saving snapshot ..."
cp *.gif "$folder"
cp *.pdf "$folder"
cp *.png "$folder"
cp *.root "$folder"
cp *.txt "$folder"
cp *.json "$folder"
cp *.tex "$folder"
cp toys/results/*.root "$folder/toys/results/"

rm *.lock

tar cjf "snapshot.tar.bz2" "$folder/"
mv "snapshot.tar.bz2" "$folder"

echo "Done."

echo "Snapshot saved in $folder"
