#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.450
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl


if [ $# -lt 1 ]; then
    echo "Need type (sim or dat 0 or 1)"
fi

TYPE=sim
if [ $1 -eq 1 ]; then
    TYPE=dat
fi

#hadd -f -j 4 "summed_${TYPE}.root" `ls output/root/*${TYPE}*fullfile*`
hadd -f -j 4 "summed_chi2file.root" `ls output/temphists/*chi2file*`
