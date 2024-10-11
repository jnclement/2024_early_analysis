#!/bin/bash
valgrind --leak-check=yes --trace-children=yes ./earlydata.sh debug 0 47289 1 1100 0
