#!/bin/bash


SUBNAME="haddmb.sub"
njob=$(( `cat chi2lists/mbfilelist.list | wc -l` / 100 ))
echo $njob
#njob=15
echo "executable = mb_hadd_condor.cmd" > $SUBNAME
#echo "concurrency_limits=CONCURRENCY_LIMIT_DEFAULT:100" >> $SUBNAME
echo "arguments = \$(Process)" >> $SUBNAME
echo "output = /sphenix/user/jocl/projects/run2024_earlydata/run/output/out/output_\$(Process)_haddmb.out" >> $SUBNAME
echo "request_cpus = 4" >> $SUBNAME
echo "should_transfer_files   = IF_NEEDED" >> $SUBNAME
echo "when_to_transfer_output = ON_EXIT" >> $SUBNAME
echo "error = /sphenix/user/jocl/projects/run2024_earlydata/run/output/err/error_\$(Process)haddmb.err" >> $SUBNAME
echo "log = /tmp/jocl_haddmb.log" >> $SUBNAME
echo "priority = 1000000000" >> $SUBNAME
echo "queue ${njob}" >> $SUBNAME

condor_submit $SUBNAME
