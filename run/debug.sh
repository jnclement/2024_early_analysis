#!/bin/bash

if [ $# -lt 5 ]; then
    echo "Need arguments nevent, run number, debug level, datorsim, segment"
    exit 1
fi

root -b -q 'run_earlydata.C("debug",'${5}','${3}','${1}','${2}',1,'${4}')'
