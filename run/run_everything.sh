#!/bin/bash

TAG=$1
NFILE=$2
RN=$3
ZS=$4
EVT=$5
C2C=$6
if [ $# -lt 6 ]; then
    echo "Need arguments (in order): tag (string), nfile (int), run number (int), ZS (int) EVT (int), C2C (int)"
    exit 1
fi

BASENAME="condor_${TAG}_${NFILE}_${RN}_${ZS}_${EVT}"
PREFIX="." #"/sphenix/user/hanpuj/test"
SUBNAME="${BASENAME}.sub"

#echo "executable = containerscripts/earlydata.sh" > $PREFIX/$SUBNAME
echo "executable = earlydata.sh" > $PREFIX/$SUBNAME
echo "arguments = ${TAG} \$(Process) ${RN} ${ZS} ${EVT} ${C2C}" >> $PREFIX/$SUBNAME
echo "priority = 100000000" >> $SUBNAME
echo "output = output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "request_memory                = 4GB" >> $PREFIX/$SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $PREFIX/$SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $PREFIX/$SUBNAME
echo "error = output/out/output_${BASENAME}_\$(Process).out" >> $PREFIX/$SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $PREFIX/$SUBNAME
echo "queue ${NFILE}" >> $PREFIX/$SUBNAME

condor_submit $PREFIX/$SUBNAME
