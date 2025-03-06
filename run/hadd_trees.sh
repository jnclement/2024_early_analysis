#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Need sample type argument"
    exit 1
fi
#rm -f output/${1}_hadd_trees/*
for i in {101..1000}; do
    echo $i
    rm lists/sim_haddlist.list
    LOWER=$(( $i * 100 ))
    UPPER=$(( $LOWER + 99 ))
    ls /sphenix/tg/tg01/jets/jocl/evt/${LOWER}/events_20250303_useE_${LOWER}_${LOWER}_${LOWER}_0.root
    if [ $? -ne 0 ]; then
	break
    fi
    for j in $(seq $LOWER $UPPER); do
	ls /sphenix/tg/tg01/jets/jocl/evt/${j}/events_20250303_useE_${1}_${j}_${j}_0.root
	if [ $? -ne 0 ]; then
	    break
	fi
	echo "/sphenix/tg/tg01/jets/jocl/evt/${j}/events_20250303_useE_${1}_${j}_${j}_0.root" >> lists/sim_haddlist.list
    done

    hadd -f -j 2 output/${1}_hadd_trees/${1}_haddtree_${i}.root `cat lists/sim_haddlist.list`
done

#hadd -n 200 -j 4 -f $1_single_tree.root `ls output/${1}_hadd_trees/*`
