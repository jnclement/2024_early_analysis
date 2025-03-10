#!/bin/bash

TAG=`echo $1 | awk -F'/' '{print $2}' | awk -F'.' '{print $1}'`
if [ $# -lt 2 ]; then
    echo "Need arguments (in order): filename (string), sampletype (string)"
    exit 1
fi

NJOB=`wc -l < $1`
NJOB=$(( ($NJOB + 9) / 10))
#NJOB=10
BASENAME="condor_${TAG}_${NJOB}_imagemaker"

SUBNAME="${BASENAME}.sub"

echo "executable = quick_jet10.sh" > $SUBNAME
echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:100" >> $SUBNAME
echo "arguments = ${1} \$(Process) ${2}" >> $SUBNAME
echo "output = output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = output/err/error_${BASENAME}_\$(Process).err" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}_\$(Process).log" >> $SUBNAME
echo "priority = 1000000000" >> $SUBNAME
echo "queue ${NJOB}" >> $SUBNAME

condor_submit $SUBNAME
