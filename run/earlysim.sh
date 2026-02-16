#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.533
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
N=5
if [ "$7" = "mb" ]; then
    N=20
fi

STARTN=$(( $2 * $N ))
NL=$(( $N - 1 ))

mkdir -p $2\_chi2
mkdir -p lists
mkdir -p ./dsts/

for i in $(seq 0 $NL); do
    SUBDIR=$(( $STARTN + $i ))
    UPLN=$(( $SUBDIR + 1 ))

    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/run_earlydata.C .
    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/${7}list/g4hits.list ./lists/g4hits.list
    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/${7}list/dst_truth_jet.list ./lists/dst_truth_jet.list
    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/${7}list/dst_calo_cluster.list ./lists/dst_calo_cluster.list
    #cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_global.list ./lists/dst_global.list
    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/${7}list/dst_mbd_epd.list ./lists/dst_mbd_epd.list
    G4HITSF=`sed -n "${UPLN}"p ./lists/g4hits.list`
    CALOCLF=`sed -n "${UPLN}"p ./lists/dst_calo_cluster.list`
    #GLOBALF=`sed -n "${UPLN}"p ./lists/dst_global.list`
    TRTHJET=`sed -n "${UPLN}"p ./lists/dst_truth_jet.list`
    DMBDEPD=`sed -n "${UPLN}"p ./lists/dst_mbd_epd.list`
    FULLTRTH=`psql FileCatalog -tA -c "select full_file_path from files where lfn = '${TRTHJET}';"`
    FULLMBEP=`psql FileCatalog -tA -c "select full_file_path from files where lfn = '${DMBDEPD}';"`
    FULLCALO=`psql FileCatalog -tA -c "select full_file_path from files where lfn = '${CALOCLF}';"`
    FULLG4HT=`psql FileCatalog -tA -c "select full_file_path from files where lfn = '${G4HITSF}';"`
    echo $CALOCLF
    #echo $GLOBALF
    #getinputfiles.pl $GLOBALF
    getinputfiles.pl $CALOCLF
    getinputfiles.pl $TRTHJET
    getinputfiles.pl $DMBDEPD
    getinputfiles.pl $G4HITSF
    #dd if="${FULLTRTH}" of="./${TRTHJET}" bs=12M
    #dd if="${FULLMBEP}" of="./${DMBDEPD}" bs=12M
    #dd if="${FULLCALO}" of="./${CALOCLF}" bs=12M
    #cp $FULLG4HT .
    #cp -r $G4HITSF ./dsts/$2/g4hits_${2}.root
    echo ""
    echo "" 
    ls
    echo ""
    echo ""
    mv $CALOCLF ./dsts/calo_cluster_${SUBDIR}.root
    #mv $GLOBALF ./dsts/$SUBDIR/global_${SUBDIR}.root
    mv $TRTHJET ./dsts/truth_jet_${SUBDIR}.root
    mv $DMBDEPD ./dsts/mbd_epd_${SUBDIR}.root
    mv $G4HITSF ./dsts/g4hits_${SUBDIR}.root
done
ls ./dsts/

ls ./dsts/calo_cluster* > ./calo_cluster_temp.list
ls ./dsts/truth_jet* > ./truth_jet_temp.list
ls ./dsts/mbd_epd* > ./mbd_epd_temp.list
ls ./dsts/g4hits* > ./g4hits_temp.list

#cp -r $TRTHJET ./dsts/$SUBDIR/truth_jet_${SUBDIR}.root
root -b -q 'run_earlydata.C("'${1}'",'${SUBDIR}',0,'${5}','${2}','${4}',0,'${6}')'
echo " "
echo " "

mkdir -p /sphenix/tg/tg01/jets/jocl/chi2/$2
for file in `ls ./$2\_chi2/`; do
    dd if="./${2}_chi2/${file}" of="/sphenix/tg/tg01/jets/jocl/chi2/${2}/${file}" bs=12M
done
