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
STARTN=$(( $2 * 20 ))
for i in {0..19}; do
    SUBDIR=$(( $STARTN + $i ))
    UPLN=$(( $SUBDIR + 1 ))
    mkdir -p $SUBDIR
    mkdir -p $2\_chi2
    mkdir -p lists
    mkdir -p ./output/smg
    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/
    mkdir -p ./dsts/$SUBDIR
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
    mv $CALOCLF ./dsts/$SUBDIR/calo_cluster_${SUBDIR}.root
    #mv $GLOBALF ./dsts/$SUBDIR/global_${SUBDIR}.root
    mv $TRTHJET ./dsts/$SUBDIR/truth_jet_${SUBDIR}.root
    mv $DMBDEPD ./dsts/$SUBDIR/mbd_epd_${SUBDIR}.root
    mv $G4HITSF ./dsts/$SUBDIR/g4hits_${SUBDIR}.root
    ls ./dsts/$SUBDIR
    #cp -r $TRTHJET ./dsts/$SUBDIR/truth_jet_${SUBDIR}.root
    root -b -q 'run_earlydata.C("'${1}'",'${SUBDIR}',0,'${5}','${2}','${4}',0,'${6}')'
    ls
    echo " "
    ls $SUBDIR
    echo " "
    ls output/smg
    mkdir -p /sphenix/tg/tg01/jets/jocl/chi2/$SUBDIR/
    for file in `ls ./$2\_chi2/`; do
	dd if="./${2}_chi2/${file}" of="/sphenix/tg/tg01/jets/jocl/chi2/${2}/${file}" bs=12M
    done
    rm -r ./$2\_chi2/*
    rm ./dsts/$SUBDIR/*
    #cp -r "./${SUBDIR}/events_${1}_${SUBDIR}_${SUBDIR}_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${SUBDIR}_${SUBDIR}_${5}.root"
done
