#!/bin/bash

for file in `ls *.imagelist`; do
    bash run_images.sh $file
done
