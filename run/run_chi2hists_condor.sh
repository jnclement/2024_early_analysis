#!/bin/bash

NJOB=`wc -l < chi2lists/chi2list.list`
#NJOB=10
BASENAME="condor_chi2hist_${NJOB}"
SUBNAME="${BASENAME}.sub"

#echo "executable = containerscripts/build_chi2hists.sh" > $SUBNAME
echo "executable = build_chi2hists.sh" > $SUBNAME
echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:200" >> $SUBNAME
echo "arguments = \$(Process)" >> $SUBNAME
echo "output = output/out/output_${BASENAME}_\$(Process).out" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = output/err/error_${BASENAME}_\$(Process).err" >> $SUBNAME
echo "log = /tmp/jocl_${BASENAME}_\$(Process).log" >> $SUBNAME
echo "priority = 10000000" >> $SUBNAME
echo "queue $NJOB" >> $SUBNAME

condor_submit $SUBNAME
