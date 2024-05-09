#!/bin/bash
#74036.1204 jocl           10/15 20:45   0+00:05:15 H  0   2198.0 dETdeta.sh 1204 20231011_nopileup 1 0 0 0

condor_q | grep jocl | while read -r line; do
    #echo $line
    HELDVAL=$(echo $line | awk '{print $6;}' )
    TAG=$(echo $line | awk '{print $11;}' )
    MCDAT=$(echo $line | awk '{print $12;}' )
    COR=$(echo $line | awk '{print $15;}' )
    NAME=$(echo $line | awk '{print $1;}' )
    NUM=$(echo $line | awk '{print $10;}' )
    #echo $HELDVAL
    if [ "$HELDVAL" = "H" ]; then
	if [ $MCDAT = 1 ]; then
	    MCDAT="mc"
	else
	    MCDAT="data"
	fi
	if [ $COR = 1 ]; then
	    COR="cor"
	else
	    COR="unc"
	fi
	FILENAME="output/evt/events_"$TAG"_"$MCDAT"_"$COR"_"$NUM".root"
	echo "removing file "$FILENAME
	rm $FILENAME
	#echo "removing job "$NAME
    fi
done
