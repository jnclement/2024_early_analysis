#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Need tag argument (string), evtnum (int), chi2check (int)"
    exit 1
fi

nmax=20000
evt=$2
c2c=$3
filecounter=0
if [ $evt -gt 100000 ]; then
    evt=0
fi
echo $evt
for rn in `ls  lists/dst_calo_run2pp*.list | awk -F'.' '{print $1}' | awk -F'/' '{print $2}' | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l lists/dst_calo_run2pp-000${rn}.list | awk '{print $1}'`
    filecounter=$(( $filecounter + $nfile ))
    if [ $filecounter -gt $nmax ]; then
	break
    fi
#    nfile=$(( ($nfile + 9) / 10 ))
    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/$rn
    mkdir -p /sphenix/tg/tg01/jets/jocl/chi2/$rn
    bash run_everything.sh $1 $nfile $rn 1 $evt $c2c
done


