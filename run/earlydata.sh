#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
export TESTINSTALL=/sphenix/user/jocl/projects/testinstall
echo $LD_LIBRARY_PATH
echo $PATH
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
SUBDIR=${3}
echo $SUBDIR
UPLN=$(( $2 + 1 ))
mkdir -p $SUBDIR
mkdir -p lists
mkdir -p /sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/
mkdir -p /sphenix/tg/tg01/jets/jocl/chi2/${SUBDIR}/
mkdir -p ./dsts/$SUBDIR
mkdir -p ./output/smg
mkdir -p output/err
mkdir -p output/out
mkdir -p $SUBDIR\_chi2
echo "Made dirs"
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/run_earlydata.C .
#cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/Fun4All_CaloDataAna.C .
echo "copied root script here"
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_calo_run2pp-000$3.list ./lists/$3.list
echo "copied dstlist here"
DSTFILE=`sed -n "${UPLN}"p "./lists/${3}.list"`
getinputfiles.pl $DSTFILE
echo "got input files"
mv $DSTFILE ./dsts/$3/${3}_${2}.root
echo "Running Fun4All now"
echo $DSTFILE
echo ./dsts/$3/${3}_${2}.root
root -b -q 'run_earlydata.C("'${1}'",'${2}',0,'${5}','${3}','${4}',1,'${6}')'
#root -b -q 'Fun4All_CaloDataAna.C(2000,"'./dsts/$3/$3_$2.root'")'
echo " "
echo " "
echo " "
echo "after run"
ls -larth
echo " "
ls -larth $SUBDIR
ls -larth $SUBDIR\_chi2/*
cp -r "./${SUBDIR}/events_${1}_${3}_${2}_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${3}_${2}_${5}.root"
cp -r ./output/smg/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/
cp -r ./output/out/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/out/
cp -r ./output/err/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/err/
cp -r ./${SUBDIR}_chi2/* /sphenix/tg/tg01/jets/jocl/chi2/$SUBDIR/
