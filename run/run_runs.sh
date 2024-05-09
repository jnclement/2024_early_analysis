#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Need tag argument (string)"
    exit 1
fi

nmax=20

for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/temp24DSTs/*0000* | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l ${rn}_ysZS.list | awk '{print $1}'`
    if [ $nfile -gt $nmax ]; then
	nfile=$nmax
    fi
    bash run_everything.sh $1 $nfile $rn 1
done

for rn in `ls  /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/tempDEST_noSZS/*0000* | awk -F'-' '{print $2}'`; do
    rn=$(expr $rn + 0)
    nfile=`wc -l ${rn}_noZS.list | awk '{print $1}'`
    if [ $nfile -gt $nmax ]; then
	nfile=$nmax
    fi
    bash run_everything.sh $1 $nfile $rn 0
done
