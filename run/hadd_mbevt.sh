#!/bin/bash

for dir in `find /sphenix/tg/tg01/jets/jocl/chi2/ -mindepth 1 -maxdepth 1 -type d`; do
    for file in `ls -U -f -1 $dir/ | grep ".*20250119.*mbtree.*"`; do
	echo $dir/$file
	echo $dir/$file >> chi2lists/mbhaddlist.list
    done
done
hadd -f -j 4 -n 100 ./output/mb_hadd/mbhadded.root `cat chi2lists/mbhaddlist.list`
