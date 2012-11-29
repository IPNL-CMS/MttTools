#! /bin/bash

echo "Sourcing GCC"
source /afs/cern.ch/sw/lcg/contrib/gcc/4.6.2/x86_64-slc5-gcc46-opt/setup.sh /afs/cern.ch/sw/lcg/contrib
gcc -v

echo "Sourcing ROOT"
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.03/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh

echo "Sourcing Python / Ganga"
export PATH=/afs/cern.ch/sw/lcg/external/Python/2.7.3/x86_64-slc5-gcc46-opt/bin/:/afs/cern.ch/sw/ganga/install/5.8.18-hotfix1/bin/:$PATH
