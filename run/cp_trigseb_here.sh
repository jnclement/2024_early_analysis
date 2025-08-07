#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"

export HOME=/sphenix/u/jocl
export TESTINSTALL=/sphenix/user/jocl/projects/testinstall
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi

SEB=$1

if [ $SEB -lt 10 ]; then
    SEB="0${1}"
fi

DSTFILE=`sed -n "${2}"p "./lists/dst_triggered_event_seb${$SEB}-*.list"`

FULLPATH=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${DSTFILE}';"`

cp $FULLPATH ./dsts/seb${1}.root
