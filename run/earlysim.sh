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
SUBDIR=$2
UPLN=$(( $2 + 1 ))
mkdir -p $SUBDIR
mkdir -p $SUBIDR\_chi2
mkdir -p lists
mkdir -p ./output/smg
mkdir -p /sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/
mkdir -p ./dsts/$SUBDIR
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/run_earlydata.C .
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/g4hits.list ./lists/g4hits.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_truth_jet.list ./lists/dst_truth_jet.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_calo_cluster.list ./lists/dst_calo_cluster.list
cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_global.list ./lists/dst_global.list
G4HITSF=`sed -n "${UPLN}"p ./lists/g4hits.list`
CALOCLF=`sed -n "${UPLN}"p ./lists/dst_calo_cluster.list`
GLOBALF=`sed -n "${UPLN}"p ./lists/dst_global.list`
TRTHJET=`sed -n "${UPLN}"p ./lists/dst_truth_jet.list`
echo $CALOCLF
echo $GLOBALF
getinputfiles.pl $GLOBALF
getinputfiles.pl $CALOCLF
#cp -r $G4HITSF ./dsts/$2/g4hits_${2}.root
echo ""
echo "" 
ls
echo ""
echo ""
mv $CALOCLF ./dsts/$2/calo_cluster_${2}.root
mv $GLOBALF ./dsts/$2/global_${2}.root
ls ./dsts/$2
#cp -r $TRTHJET ./dsts/$2/truth_jet_${2}.root
root -b -q 'run_earlydata.C("'${1}'",'${2}',0,'${5}','${2}','${4}',0,'${6}')'
ls
echo " "
ls $SUBDIR
echo " "
ls output/smg
cp -r "./${SUBDIR}/events_${1}_${2}_${2}_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${2}_${2}_${5}.root"
cp ./output/smg/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/
