#!/bin/bash


source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
hadd "output/sumroot/summed_${1}_${2}_base.root" `ls output/root/*${1}*${2}*fullfile*`
