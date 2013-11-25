#! /bin/bash

echo "Sourcing root"
source /gridgroup/cms/brochet/root/root-v5-34/bin/thisroot.sh

echo "Sourcing gcc / gdb / python / boost"
source /gridsoft/cmss/slc5_amd64_gcc472/external/gcc/4.7.2/etc/profile.d/init.sh
source /gridsoft/cmss/slc5_amd64_gcc472/external/gdb/7.5-cms2/etc/profile.d/init.sh
source /gridsoft/cmss/slc5_amd64_gcc472/external/python/2.7.3-cms6/etc/profile.d/init.sh
source /gridsoft/cmss/slc5_amd64_gcc472/external/boost/1.51.0-cms4/etc/profile.d/init.sh

export BOOST_ROOT
