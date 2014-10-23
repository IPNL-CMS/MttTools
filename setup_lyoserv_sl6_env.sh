#! /bin/bash

echo "Sourcing root"
#source /gridgroup/cms/brochet/root/root-combine-sl6/bin/thisroot.sh
source /gridgroup/cms/brochet/root/root-sl6/bin/thisroot.sh

echo "Sourcing gcc / gdb / python / boost"
source /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gcc/4.7.2-cms/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/gdb/7.6-cms/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/python/2.6.4/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/boost/1.51.0-cms5/etc/profile.d/init.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/git/1.8.3.1-cms2/etc/profile.d/init.sh
export GIT_EXEC_PATH=/cvmfs/cms.cern.ch/slc6_amd64_gcc472/external/git/1.8.3.1-cms2/libexec/git-core

export BOOST_ROOT
