#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Need tag argument (string), evtnum (int), chi2check (int)"
    exit 1
fi
rn=0
evt=$2
c2c=$3
nfile=10000
for i in {0..10000}; do
    mkdir -p /sphenix/tg/tg01/jets/jocl/evt/$i
    mkdir -p /sphenix/tg/tg01/jets/jocl/err/$i
    mkdir -p /sphenix/tg/tg01/jets/jocl/out/$i
done
if [ $evt -gt 1000 ]; then
    evt=0
fi
bash run_everysim.sh $1 $nfile $rn 1 $evt $c2c

#for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/tempDEST_noSZS/*0000* | awk -F'-' '{print $2}'`; do
#    rn=$(expr $rn + 0)
#    nfile=`wc -l ${rn}_noZS.list | awk '{print $1}'`
#    if [ $nfile -gt $nmax ]; then
#	nfile=$nmax
#    fi
#    bash run_everything.sh $1 $nfile $rn 0
#done
