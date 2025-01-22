#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.450
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
mkdir -p lists
mkdir -p output/simhists
cp /sphenix/user/jocl/projects/run2024_earlydata/run/lists/sim.imagelist .

NJOB=$(( ($2 * 10) + 1))
mkdir -p treefiles_$NJOB
FILENAME=lists/$NJOB.imagelist
cat sim.imagelist | tail -n +$NJOB | head -n 10 > ./output/templist_$NJOB.list

cat ./output/templist_$NJOB.list

while read -r evtfile; do
    echo $evtfile
    cp $evtfile ./treefiles_$NJOB
done < ./output/templist_$NJOB.list

ls treefiles_$NJOB/* > $FILENAME
cp /sphenix/user/jocl/projects/run2024_earlydata/run/quick_jet10.C .
cp /sphenix/user/jocl/projects/run2024_earlydata/run/dlUtility.h .
root -b -q 'quick_jet10.C("'${FILENAME}'",'${2}',1)'
cp -r output/simhists/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/simhists/
