#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi
STARTN=$(( $2 * 10 ))
for i in {0..9}; do
    SUBDIR=$(( $STARTN + $i ))
    UPLN=$(( $SUBDIR + 1 ))
    mkdir -p $SUBDIR
    mkdir -p $2\_chi2
    mkdir -p lists
    mkdir -p ./output/smg
    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/
    mkdir -p ./dsts/$SUBDIR
    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/run_cosmicsim.C .
    cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_cosmics.list ./lists/dst_cosmics.list
    CALOCLF=`sed -n "${UPLN}"p ./lists/dst_cosmics.list`

    FULLCALO=$CALOCLF
    echo $CALOCLF
    echo ""
    echo "" 
    ls
    echo ""
    echo ""
    cp $FULLCALO ./dsts/$SUBDIR/cosmics_${SUBDIR}.root
    ls ./dsts/$SUBDIR
    #cp -r $TRTHJET ./dsts/$SUBDIR/truth_jet_${SUBDIR}.root
    root -b -q 'run_cosmicsim.C("'${1}'",'${SUBDIR}',0,'${3}','${2}',0,0,1)'
    ls
    echo " "
    ls $SUBDIR
    echo " "
    ls output/smg
    mkdir -p /sphenix/tg/tg01/jets/jocl/chi2/$SUBDIR/
    cp -r ./$2\_chi2/* /sphenix/tg/tg01/jets/jocl/chi2/$SUBDIR/
    rm -r ./$2\_chi2/*
    #cp -r "./${SUBDIR}/events_${1}_${SUBDIR}_${SUBDIR}_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${SUBDIR}_${SUBDIR}_${5}.root"
done
