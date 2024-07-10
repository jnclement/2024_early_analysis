#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Need tag argument (string), evtnum (int), chi2check (int)"
    exit 1
fi

nmax=1500
evt=$2
c2c=$3
if [ $evt -gt 10000 ]; then
    evt=0
fi
echo $evt
for rn in `ls  *.list | awk -F'.' '{print $1}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l ${rn}.list | awk '{print $1}'`
    if [ $nfile -gt $nmax ]; then
	nfile=$nmax
    fi
#    nfile=$(( ($nfile + 9) / 10 ))
    bash run_everything.sh $1 $nfile $rn 1 $evt $c2c
done

#for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/tempDEST_noSZS/*0000* | awk -F'-' '{print $2}'`; do
#    rn=$(expr $rn + 0)
#    nfile=`wc -l ${rn}_noZS.list | awk '{print $1}'`
#    if [ $nfile -gt $nmax ]; then
#	nfile=$nmax
#    fi
#    bash run_everything.sh $1 $nfile $rn 0
#done
