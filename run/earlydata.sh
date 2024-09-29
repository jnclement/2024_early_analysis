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
UPLN=$(( $2 + 1 ))
mkdir -p $SUBDIR
mkdir -p lists
mkdir -p ./output/smg
mkdir -p /sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/
mkdir -p ./dsts/$SUBDIR

cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/run_earlydata.C .
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/$3.list ./lists/$3.list
DSTFILE=`sed -n "${UPLN}"p "./lists/${3}.list"`
cp -r $DSTFILE ./dsts/$3/${3}_${2}.dst

root -b -q 'run_earlydata.C("'${1}'",'${2}',0,'${5}','${3}','${4}',1,'${6}')'
ls
echo " "
ls $SUBDIR
echo " "
ls output/smg
cp -r "./${SUBDIR}/events_${1}_${3}_${2}_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${3}_${2}_${5}.root"
cp ./output/smg/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/
