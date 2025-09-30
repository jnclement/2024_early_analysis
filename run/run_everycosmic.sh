#!/bin/bash

TAG=$1
NFILE=$2
EVT=$3
if [ $# -lt 3 ]; then
    echo "Need arguments (in order): tag (string), nfile (int), EVT (int)"
    exit 1
fi

BASENAME="condor_${TAG}_${NFILE}_${EVT}"

SUBNAME="${BASENAME}.sub"

echo "executable = cosmicsim.sh" > $SUBNAME
echo "arguments = ${TAG} \$(Process) ${EVT}" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/run2024_earlydata/run/output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "request_memory                = 4GB" >> $SUBNAME
echo "error = /sphenix/user/jocl/projects/run2024_earlydata/run/output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $SUBNAME
echo "queue ${NFILE}" >> $SUBNAME

condor_submit $SUBNAME
