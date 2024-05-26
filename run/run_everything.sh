#!/bin/bash

TAG=$1
NFILE=$2
RN=$3
ZS=$4
EVT=$5
if [ $# -lt 5 ]; then
    echo "Need arguments (in order): tag (string), nfile (int), run number (int), ZS (int) EVT (int)"
    exit 1
fi

BASENAME="condor_${TAG}_${NFILE}_${RN}_${ZS}_${EVT}"

SUBNAME="${BASENAME}.sub"

echo "executable = earlydata.sh" > $SUBNAME
echo "arguments = ${TAG} \$(Process) ${RN} ${ZS} ${EVT}" >> $SUBNAME
echo "output = output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = output/err/error_${BASENAME}_\$(Process).err" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $SUBNAME
echo "queue ${NFILE}" >> $SUBNAME

condor_submit $SUBNAME
