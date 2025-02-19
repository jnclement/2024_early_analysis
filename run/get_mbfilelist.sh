#!/bin/bash

mkdir -p chi2lists
rm chi2lists/mbfilelist.list
for dir in `find /sphenix/tg/tg01/jets/jocl/chi2/ -mindepth 1 -maxdepth 1 -type d`; do
    echo $dir
    for file in `ls $dir -U -f -1 | grep mbtree | grep 20250120`; do
	name=`echo $dir/$file`
	echo $name >> chi2lists/mbfilelist.list
    done
done
