#!/bin/bash

TAG=$1
if [ $# -lt 1 ]; then
    echo "Need arguments (in order): filename"
    exit 1
fi

NJOB=`wc -l < $TAG`
NJOB=$(( ($NJOB + 9) / 10))
BASENAME="condor_${TAG}_${NJOB}_imagemaker"

SUBNAME="${BASENAME}.sub"

echo "executable = quickroot.sh" > $SUBNAME
echo "arguments = ${TAG} \$(Process)" >> $SUBNAME
echo "output = /sphenix/tg/tg01/jets/jocl/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = /sphenix/tg/tg01/jets/jocl/err/error_${BASENAME}_\$(Process).err" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}_\$(Process).log" >> $SUBNAME
echo "queue ${NJOB}" >> $SUBNAME

condor_submit $SUBNAME
