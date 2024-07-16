#!/bin/bash

hadd "output/sumroot/summed_${1}_${2}_base.root" `ls output/root/*${1}*${2}*fullfile*`
