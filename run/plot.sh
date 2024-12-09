#!/bin/bash

#source /opt/sphenix/core/bin/sphenix_setup.sh -n
#source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
#export HOME=/sphenix/u/jocl

root -b -q 'plot.C("summed_dat'$1'.root","summed_sim.root","lists/sumdatlist'$1'.list","rmg'$1'")' #,'${1}')'
