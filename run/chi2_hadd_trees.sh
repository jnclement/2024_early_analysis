#!/bin/bash

for rn in `cat fullgoodrunlist.list`; do
    hadd -j 2 -f output/chi2_hadd_trees/$rn\_haddtree.root `ls /sphenix/tg/tg01/jets/jocl/chi2/${rn}/*20250303_useE_forHanpu*`
