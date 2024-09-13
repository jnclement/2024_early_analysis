#!/bin/bash

for file in `ls lists/*.imagelist`; do
    for file2 in `cat $file`; do
	find $file2 -mtime +2 -exec rm {} \;
    done
done
