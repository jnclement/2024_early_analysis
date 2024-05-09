#!/bin/bash

for basefile in `basename -a /home/jocl/datatemp/*thirdgo*`; do
    root -b -q 'quickroot.C("'$basefile'")'
done
