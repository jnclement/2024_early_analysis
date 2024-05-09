#!/bin/bash

if [ $2 -eq 1 ]; then
    ls /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/temp24DSTs/*$1* > $1_ysZS.list
fi

if [ $2 -eq 0 ]; then
    ls /sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/tempDEST_noSZS/*$1* > $1_noZS.list
fi
#ls /sphenix/tg/tg01/commissioning/CaloCalibWG/dlis/DST_PRDF-00041725-0000.root > list.list
