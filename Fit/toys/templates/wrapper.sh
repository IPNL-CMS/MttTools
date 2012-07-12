#! /bin/sh

# ROOT environment
. /gridsoft/ipnls/env/root_new64.sh

# GCC
source /gridsoft/cmss/slc5_amd64_gcc462/external/gcc/4.6.2/etc/profile.d/init.sh

# cd to right directory
cd $4

# define variables
mass=$1
num_exp=$2
job_index=$3
btag=$7

# launch analysis
./fitMtt -m $mass --limit-curve --toys ${num_exp} --index ${job_index} --path \"$5\" --output-path \"$6\" --no-text-files --no-figs --no-root-files --b-tag $7 --batch

error=$?

if [ $error != 0 ]; then
  exit $error
fi

chmod a+rw $6/*job$3.root

exit $error
