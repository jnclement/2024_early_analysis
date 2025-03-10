#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
NJOB=$(( $1 * 10 ))
NJOB=$(( $NJOB + 1 ))
NJUP=$(( $NJOB + 10 ))
mkdir -p lists
mkdir -p put
mkdir -p output/err
mkdir -p output/out
for i in $(seq $NJOB $NJUP); do
    name=`tail -n +$i /sphenix/user/jocl/projects/run2024_earlydata/run/chi2lists/chi2list.list | head -n 1`
    echo $name
    cp $name.root ./put
    bn=`basename $name`
    name2=`ls put/$bn.root | awk -F'.' '{print $1}'`
    echo $name2
    rn=`ls $name2.root | awk -F'_' '{print $5}'`
    cp /sphenix/user/jocl/projects/run2024_earlydata/run/build_chi2hists.C .
    cp /sphenix/user/jocl/projects/run2024_earlydata/run/dlUtility.h .
    root -b -q 'build_chi2hists.C("'$name2'",'$rn')'
    cp -r $name2\_hist.root /sphenix/user/jocl/projects/run2024_earlydata/run/output/temphists/
done
