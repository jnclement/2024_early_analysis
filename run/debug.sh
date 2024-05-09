#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Need arguments nevent, run number"
    exit 1
fi

root -b -q 'run_earlydata.C("debug",0,1,'${1}','${2}',1)'
