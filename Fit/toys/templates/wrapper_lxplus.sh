#! /bin/sh

echo "Sourcing GCC"
source /afs/cern.ch/sw/lcg/contrib/gcc/4.6.2/x86_64-slc5-gcc46-opt/setup.sh /afs/cern.ch/sw/lcg/contrib
gcc -v

echo "Sourcing ROOT"
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.03/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh

# cd to right directory
cd $5

# define variables
mass=$2
num_exp=$3
job_index=$4
btag=$8

# launch analysis
./fitMtt -i $1 -m $mass --limit-curve --toys ${num_exp} --index ${job_index} --path \"$6\" --output-path \"$7\" --no-text-files --no-figs --no-root-files --b-tag $8 --batch

error=$?

if [ $error != 0 ]; then
  exit $error
fi

chmod a+rw $6/*job$3.root

exit $error
