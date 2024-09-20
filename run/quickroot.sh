#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
mkdir -p lists
mkdir -p output/root

IMN=`echo $1 | awk -F'.' '{print $1}' | awk -F'/' '{print $2}'`

NJOB=$(( ($2 * 10) + 1))
mkdir -p treefiles_$IMN\_$NJOB
FILENAME=lists/$IMN\_$NJOB.imagelist
cp `ls /sphenix/tg/tg01/jets/jocl/evt/$IMN/* | tail -n +$NJOB | head -n 10` treefiles_$IMN\_$NJOB

ls treefiles_$IMN\_$NJOB/* > $FILENAME
cp /sphenix/user/jocl/projects/run2024_earlydata/run/quickroot.C .
cp /sphenix/user/jocl/projects/run2024_earlydata/run/dlUtility.h .
root -b -q 'quickroot.C("'${FILENAME}'",0)'
cp -r output/root/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/root
