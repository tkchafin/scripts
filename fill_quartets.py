#!/usr/bin/python

import os
import sys
import itertools
import collections

if len(sys.argv) < 3:
	print("Usage: fill_quartets.py <CF_file> <taxon list>")

CF=sys.argv[1]
all=sys.argv[2]

spoof=True #hard coded option

#list of lists, capturing sampled quartets
sampled = list()

with open(CF, 'r') as fh:
	try:
		seen = list()
		for line in fh:
			if not line:
				continue
			else:
				stuff = line.split(",")
				seen = sorted(stuff[0:4])
				sampled.append(seen)
	except IOError:
		print("Could not read file ",CF)
		sys.exit(1)
	finally:
		fh.close()

all_quartets=list()
all_tax = list()
with open(all, 'r') as fh:
	try:
		all = list()
		for line in fh:
			line=line.strip()
			if not line:
				continue
			else:
				all_tax.append(line)
	except IOError:
		print("Could not read file ",CF)
		sys.exit(1)
	finally:
		fh.close()

all_comb = list(itertools.combinations(all_tax,4))
for comb in all_comb:
	all_quartets.append(sorted(list(comb)))

#print("Writing all missing quartets to stdout...")

for quartet in all_quartets:
	miss=True
	for sample in sampled:
		if set(quartet) == set(sample):
			miss=False
	if miss==True:
		if spoof:
			oline = "";
			for tax in quartet:
				oline = oline + str(tax) + ","
			oline = oline + "0.333333333333334,0.333333333333333,0.333333333333333"
			print(oline)
		else:
			print(quartet)
