#!/bin/bash

pre=$1

mkdir $pre-output

cd $pre-output

echo "Sumchains...."
pyco-sumchains -s 100 ../$pre-state-run-*.log &> $pre-sumchains.txt

echo "Getting optimal number for burnin..."
ch=$pre"-sumchains.txt"
samps=`grep "maximized" $ch | sed 's/.*: //g' | sed 's/ samples.*//g'`

echo "Removing $samps samples!"

echo "sumcoevolity..."
yam=$pre".yaml"
p=$pre"-"
sumcoevolity -b $samps -c ../$yam -p $p -n 1000000 ../$pre-state-run*.log

echo "pyco-sumevents...."
pyco-sumevents -p $p -f $pre-sumcoevolity-results-nevents.txt

echo "pyco-sumtimes..."
pyco-sumtimes -p $p -f -b $samps -z ../$pre-state*.log

