#!/bin/bash

if [ $# -lt 4 ]; then
    echo "Need tag argument (string), evtnum (int), chi2check (int), dowf (int)"
    exit 1
fi

nmax=10000
evt=$2
c2c=$3
dowf=$4
filecounter=0
if [ $evt -gt 100000 ]; then
    evt=0
fi
echo $evt
for rn in `cat listrunnumber.txt`; do #`cat newgoodrunlist.list`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l lists/dst_jetcalo-000${rn}.list | awk '{print $1}'`
    njob=$(( $nfile + 63 ))
    njob=$(( $njob / 64 ))
    filecounter=$(( $filecounter + $njob ))
    
#    nfile=$(( ($nfile + 9) / 10 ))
    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/$rn
    mkdir -p /sphenix/tg/tg01/jets/jocl/chi2/$rn
    echo $rn $filecounter
    echo $evt
    bash run_everything.sh $1 $njob $rn 1 $evt $c2c $dowf
    if [ $filecounter -gt $nmax ]; then
	break
    fi
done


