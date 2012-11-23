#!/bin/bash

working_dir=`pwd`/works/
input_dir=`pwd`/inputs/
output_dir=$1

analysisName=`cat ../analysis.json | grep "\"name\"" | awk -F' ' '{print $2;}' | tr -d \",`

function test_root_file()
{
  file=$1
  code=$(cat <<EOF
  {
    TFile *f = TFile::Open("${file}");
    if (! f) {
      return 1;
    } else {
    f->Close();
    delete f;
    return 0;
  }
}
EOF
)

echo "$code" > root_file_test.C
root -b -q -n -l root_file_test.C &> /dev/null
error=$?
[ -f root_file_test.C ] && rm root_file_test.C
return $error
}

echo "Merging..."
for mass in $(cat liste_mass.txt); do
  filename="${output_dir}data_2012_nominal_${mass}_toylimit_${analysisName}_*.root"
  merged_filename="${output_dir}data_2012_nominal_${mass}_toylimit_${analysisName}.root"
  hadd "${merged_filename}" ${filename}
  error=$?
  if [ $error == 0 ]; then
    if [ -f "${merged_filename}" ]; then
      test_root_file "${merged_filename}"
      error=$?
      if [ $error == 0 ]; then
        # We can safely remove files
        echo "Merged file is a valid ROOT file. Removing sub-files..."
        rm -f ${filename}
        echo "Done."
      fi
    fi
  fi
done
echo "Done!"
