#! /bin/sh

# ROOT environment
. /gridsoft/ipnls/env/root_new64.sh

# GCC
source /gridsoft/cmss/slc5_amd64_gcc462/external/gcc/4.6.2/etc/profile.d/init.sh

# GDB
source /gridsoft/cmss/slc5_amd64_gcc462/external/gdb/7.3.1/etc/profile.d/init.sh

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

chmod a+rw $7/*job$3.root

exit $error
