#!/bin/bash

#if [ $2 -eq 1 ]; then
#    ls /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/temp24DSTs/*$1* > $1_ysZS.list
#fi



#if [ $2 -eq 0 ]; then
#    ls /sphenix/lustre01/sphnxpro/commissioning/slurp/caloy2test/run_00042900_00043000/*$1* > $1.list
#fi
#ls /sphenix/tg/tg01/commissioning/CaloCalibWG/dlis/DST_PRDF-00041725-0000.root > list.list

UPPER=`awk -v n=$1 'BEGIN{print int((n+50)/100) * 100}'`
if [ $UPPER -lt $1 ]; then
    UPPER=$(($UPPER + 100))
fi

LOWER=$((UPPER - 100))

if [ $LOWER -lt 100000 ]; then
    LOWER=000$LOWER
else
    LOWER=00$LOWER
fi

if [ $UPPER -lt 100000 ]; then
    UPPER=000$UPPER
else
    UPPER=00$UPPER
fi

dir="physics"
echo $UPPER
if [ $2 -eq 0 ]; then
    dir="commissioning"
fi

ls /sphenix/lustre01/sphnxpro/$dir/slurp/caloy2test/run_$LOWER\_$UPPER/*$1* > $1.list
