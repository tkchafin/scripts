#!/bin/bash 

#Tyler K. Chafin 
#July 23 2021
#Generates n bootstrap samples of an input newick-formatted file of trees
#Email: tylerkchafin@gmail.com with issues

if [ $1 ] && [ $2 ]; 
then 
  trees="$1"
  n=$2
else 
  printf "\nUsage: $0 <tree file> <n>\n\n"
  exit 1
fi

for i in `seq 1 $n`; 
do
  ofile="b_"$i".tre"
  shuf -r -n $n $trees > $ofile
done