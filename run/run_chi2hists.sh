#!/bin/bash

mkdir -p chi2lists
rm chi2lists/chi2list.list
for file in `ls -U -f -1 /sphenix/tg/tg01/jets/jocl/chi2/*newbigprod*`; do
    name=`echo $file | awk -F"." '{print $1}'`
    echo $name >> chi2lists/chi2list.list
done
bash run_chi2hists_condor.sh
