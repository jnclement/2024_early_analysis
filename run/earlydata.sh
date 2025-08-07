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
STARTN=$(( $2 * 5 ))
echo $SUBDIR

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
#cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_jetcalo_run2pp-000$3.list ./lists/$3\_jetcalo.list
echo "copied dstlist here"

cp -r /sphenix/user/jocl/projects/run2024_earlydata/run/lists/dst_calofitting_run2pp-000$3.list ./lists/$3.list
echo "copied dstlist here"

for i in {0..4}; do
    UPLN=$(( $STARTN + $i + 1 ))
    #DSTFILE=`sed -n "${UPLN}"p "./lists/${3}_jetcalo.list"`
    #if [ -z "${DSTFILE}" ]; then
    #exit 0
    #fi
    #getinputfiles.pl $DSTFILE
    #FULLPATH=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${DSTFILE}';"`
    #cp $FULLPATH ./$DSTFILE
    #echo "got input files"
    #mv $DSTFILE ./dsts/$3/${3}_$UPLN\_jetcalo.root

    DSTFILE=`sed -n "${UPLN}"p "./lists/${3}.list"`
    if [ -z "${DSTFILE}" ]; then
	exit 0
    fi
    #getinputfiles.pl $DSTFILE
    FULLPATH=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${DSTFILE}';"`
    cp $FULLPATH ./$DSTFILE
    echo "got input files"
    mv $DSTFILE ./dsts/$3/${3}_$UPLN.root
    echo "Running Fun4All now"
    echo $DSTFILE
    echo ./dsts/$3/${3}_$UPLN.root
    root -b -q 'run_earlydata.C("'${1}'",'$UPLN',0,'${5}','${3}','${4}',1,'${6}')'
    echo " "
    echo " "
    echo " "
    echo "after run"
    ls -larth
    echo " "
    ls -larth $SUBDIR
    ls -larth $SUBDIR\_chi2/*
    #cp -r "./${SUBDIR}/events_${1}_${3}_$UPLN_${5}.root" "/sphenix/tg/tg01/jets/jocl/evt/${SUBDIR}/events_${1}_${3}_$UPLN_${5}.root"
    rm ./dsts/$3/${3}_$UPLN.root
    cp -r ./output/smg/* /sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/
    cp -r ./${SUBDIR}_chi2/* /sphenix/tg/tg01/jets/jocl/chi2/$SUBDIR/
    rm ./output/smg/*
    rm ./${SUBDIR}_chi2/*
done

