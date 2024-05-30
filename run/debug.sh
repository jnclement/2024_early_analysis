#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Need arguments nevent, run number, debug level, datorsim"
    exit 1
fi

root -b -q 'run_earlydata.C("debug",0,'${3}','${1}','${2}',1,'${4}')'
