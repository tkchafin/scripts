#!/bin/bash

if [[ $# -eq 0 ]]; then 
	echo "Usage: sumls [-m units in MB][-g units in GB][-k units in KB][/path/to/ls]" 
	exit 0;
fi;

case $1 in 
	-m) mult=.000000953674; unit=MB ;;
	-g) mult=.00000000093192; unit=GB ;; 
	-k) mult=.000976563; unit=KB ;;
esac; 
 

pre_num=`ls -lR $2 | awk '{sum+=$5}END{print sum}'`;
adj_num=`echo "($pre_num*$mult)/1; scale=3" | bc`;
echo "There are $adj_num $unit in $0"; 


