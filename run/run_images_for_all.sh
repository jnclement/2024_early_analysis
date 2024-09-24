#!/bin/bash

for file in `ls lists/4*.imagelist`; do
    bash run_images.sh $file
done
