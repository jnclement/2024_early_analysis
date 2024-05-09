#!/bin/bash

TAG=$1
if [ $# -lt 1 ]; then
    echo "Need arguments (in order): filename"
    exit 1
fi

BASENAME="condor_${TAG}_imagemaker"

SUBNAME="${BASENAME}.sub"

echo "executable = quickroot.sh" > $SUBNAME
echo "arguments = ${TAG}" >> $SUBNAME
echo "output = output/out/output_${BASENAME}.out" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = output/err/error_${BASENAME}.err" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}.log" >> $SUBNAME
echo "queue 1" >> $SUBNAME

condor_submit $SUBNAME
