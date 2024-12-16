#!/bin/bash

mkdir -p chi2lists
rm chi2lists/chi2list.list
for dir in `find /sphenix/tg/tg01/jets/jocl/chi2/ -mindepth 1 -maxdepth 1 -type d`; do
    for file in `ls -U -f -1 $dir/* | grep ".*goodcuts.*chi2file.*"`; do
	name=`echo $file | awk -F"." '{print $1}'`
	echo $name >> chi2lists/chi2list.list
    done
done
bash run_chi2hists_condor.sh
