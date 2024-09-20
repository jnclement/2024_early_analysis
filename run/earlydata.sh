#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
SUBDIR=${3}
mkdir -p $SUBDIR
mkdir -p lists
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/run_earlydata.C .
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/$3.list ./lists/$3.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/g4hits.list ./lists/g4hits.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_truth_jet.list ./lists/dst_truth_jet.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_calo_waveform.list ./lists/dst_calo_waveform.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_global.list ./lists/dst_global.list
root -b -q 'run_earlydata.C("'${1}'",'${2}',0,'${5}','${3}','${4}',1,'${6}')'
ls
echo " "
ls $SUBDIR
cp -r "./${SUBDIR}/events_${1}_${3}_${2}_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${3}_${2}_${5}.root"
