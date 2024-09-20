#!/bin/bash


SUBNAME="hadd.sub"

echo "executable = haddfinal.sh" > $SUBNAME
echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:12500" >> $SUBNAME
echo "arguments = ${1}" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/run2024_earlydata/run/output/out/output_hadd.out" >> $SUBNAME
echo "request_cpus = 4" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = /sphenix/user/jocl/projects/run2024_earlydata/run/output/err/error_hadd.err" >> $SUBNAME
echo "log = /tmp/jocl_hadd.log" >> $SUBNAME
echo "queue 1" >> $SUBNAME

condor_submit $SUBNAME
