#!/bin/bash

CreateDstList.pl DST_CALO_run2pp --build ana450 --cdb 2024p009 --run $1
mv dst_calo_run2pp-000$1.list lists
