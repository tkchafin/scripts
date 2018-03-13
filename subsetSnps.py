#!/usr/bin/python

import sys

print("Usage: ",sys.argv[0], " SubsetFile StructureFile")

print(sys.argv[1])

#Get loci to keep
loci = []
with open(sys.argv[1]) as file_object:
	count = 0
	for line in file_object:
		if count == 0:
			count += 1
			continue
		line = line.strip()
		t = line.split()
		loc = t[0].replace("\"","")
		loci.append(int(loc))
file_object.close()

output = open("out.str", "w")
with open(sys.argv[2]) as file_2:
	for line in file_2:
		line = line.strip()
		t = line.split()
		#print("Line: ", t[0])
		col = 0
		snp = 0
		for c in t:
			if col in (0,1):
				#print(c)
				if col == 0:
					output.write(c)
				else:
					stuff = "\t" + c
					output.write(stuff)
				col += 1
			else:
				snp += 1
				if snp in loci:
					#print(snp)
					stuff = "\t" + c
					output.write(stuff)
		output.write("\n")
output.close()
