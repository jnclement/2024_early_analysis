#!/bin/bash


source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl

hadd -j 4 -f /sphenix/user/jocl/projects/run2024_earlydata/run/output/chi2_haddhists/chi2_hadded_$1.root `cat /sphenix/user/jocl/projects/run2024_earlydata/run/chi2_hadd_condor_files.txt | head -n $(( $(( $1 + 1 )) * 100 )) | tail -n 100`
