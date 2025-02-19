#!/bin/bash


SUBNAME="hadd.sub"
find /sphenix/user/jocl/projects/run2024_earlydata/run/output/simhists/ -type f -name '*run_mb_dat*.root' > mbmc_hadd_condor_files.txt
njob=$(( `cat mbmc_hadd_condor_files.txt | wc -l` / 100 ))
#njob=15
echo "executable = mbmc_hadd_condor.cmd" > $SUBNAME
#echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:100" >> $SUBNAME
echo "arguments = \$(Process)" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/run2024_earlydata/run/output/out/output_\$(Process)_hadd.out" >> $SUBNAME
echo "request_cpus = 4" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = /sphenix/user/jocl/projects/run2024_earlydata/run/output/err/error_\$(Process)hadd.err" >> $SUBNAME
echo "log = /tmp/jocl_hadd.log" >> $SUBNAME
echo "queue ${njob}" >> $SUBNAME

condor_submit $SUBNAME
