#!/bin/bash

if [ $# -lt 1 ]; then
    echo "need tag! usually the date as yyyymmdd."
    exit 1
fi

bash run_runs_sim.sh $1_jet60 0 1 jet60
bash run_runs_sim.sh $1_jet50 0 1 jet50
bash run_runs_sim.sh $1_jet40 0 1 jet40
bash run_runs_sim.sh $1_jet30 0 1 jet30
bash run_runs_sim.sh $1_jet20 0 1 jet20
bash run_runs_sim.sh $1_jet12 0 1 jet15
#bash run_runs_sim.sh $1_jet10 0 1 jet10
bash run_runs_sim.sh $1_jet5 0 1 jet5
#bash run_runs_sim.sh $1_mb 0 1 mb
