#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Need tag argument (string), evtnum (int), chi2check (int)"
    exit 1
fi

nmax=350
evt=$2
c2c=$3
if [ $evt -gt 10000 ]; then
    evt=0
fi
echo $evt
for rn in `ls  lists/*.list | awk -F'.' '{print $1}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l lists/${rn}.list | awk '{print $1}'`
    if [ $nfile -gt $nmax ]; then
	nfile=$nmax
    fi
#    nfile=$(( ($nfile + 9) / 10 ))
    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/$rn
    mkdir -p /sphenix/tg/tg01/jets/jocl/err/$rn
    mkdir -p /sphenix/tg/tg01/jets/jocl/out/$rn
    bash run_everything.sh $1 $nfile $rn 1 $evt $c2c
done


