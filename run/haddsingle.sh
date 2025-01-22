#!/bin/bash


source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.450
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
hadd -f "output/sumroot3/summed_${1}_${2}_base.root" `ls output/root3/*${1}*${2}*fullfile*`
