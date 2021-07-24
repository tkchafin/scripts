#!/bin/bash 

#Tyler K. Chafin 
#July 23 2021
#Generates n bootstrap samples of an input newick-formatted file of trees
#Email: tylerkchafin@gmail.com with issues

if [ $1 ]; 
then 
  fasta="$1"
else 
  printf "\nUsage: $0 <fasta file>\n\n"
  exit 1
fi

#calculate median sequence length
median=`grep -v ">" $fasta | awk 'BEGIN{FS=""}{print NF}' | sort -n | awk '{a[NR]=$0}END{print(NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2}'`

#select out sequences equal to or above median length 
grep -B1 "^[A-Za-z]\{$median,\}" $fasta | sed "/^--$/d" > $fasta".filter"
