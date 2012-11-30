#! /bin/bash

echo "Sourcing root"
source /gridsoft/ipnls/env/root_new64.sh

echo "Sourcing gcc / gdb"
source /gridsoft/cmss/slc5_amd64_gcc462/external/gcc/4.6.2/etc/profile.d/init.sh;
source /gridsoft/cmss/slc5_amd64_gcc462/external/gdb/7.3.1/etc/profile.d/init.sh;
source /gridsoft/cmss/slc5_amd64_gcc462/external/python/2.6.4/etc/profile.d/init.sh

echo "Sourcing python"
export PATH=/gridsoft/ipnls/python/Python-2.7.2/:$PATH
