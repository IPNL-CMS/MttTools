./process_all --data --b-tag 2
./process_all --semimu --b-tag 2
./process_all --semie --b-tag 2
mv *Histo*.root 2btag
./process_all --data --b-tag 1
./process_all --semimu --b-tag 1
./process_all --semie --b-tag 1
mv *Histo*.root 1btag
