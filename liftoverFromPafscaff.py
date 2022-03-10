#!/usr/bin/python

import sys
import os
import getopt
import csv
import pandas as pd
import functools

#sorry this code isn't really commented
#i'm in a hurry
#and tired

def main():
	params = parseArgs()
	mapper=dict()
	coords=pd.read_csv(params.coords, sep="\t", names=["scaffold", "scaffold_pos"])

	#capture mappings from pafscaff headers
	with open(params.paf, "r") as fh:
		for line in fh:
			l=line.split()
			chr=l[0].split(".")[0].replace(">","")
			if l[1] == "RevComp":
				revcomp=True
				scaffold=l[2]
				start=int(l[7].split(":")[0].replace(",",""))
				end=int(l[7].split(":")[1].replace(";","").replace(",",""))
			else:
				revcomp=False
				scaffold=l[1]
				#print(l)
				#print(line[6])
				start=int(l[6].split(":")[0].replace(",",""))
				end=int(l[6].split(":")[1].replace(";","").replace(",",""))
			mapper[scaffold]=[chr, revcomp, start, end]
			#print(scaffold, ":", mapper[scaffold])

			#sys.exit()

	#...watch out for off-by-one errors
	def liftover(mapper, row):
		#print(row)
		if row[0] in mapper:
			convert = mapper[row[0]]
			#print(convert)
			if convert[1]:
				#revcomp
				new_coord=convert[3]-row[1]
			else:
				#not revcomp
				new_coord=convert[2]+row[1]
			row['chrom']=convert[0]
			row['chrom_pos']=new_coord
		else:
			#print("Scaffold",str(row[0]), "not placed in pafscaff output"
			row['chrom'] = "NA"
			row['chrom_pos'] = 0
		return(row)
			#return(["NA", 0])
	#print(coords)
	liftover_call = functools.partial(liftover, mapper)
	coords=coords.apply(liftover_call, axis = 1)
	print(coords)

	coords.to_csv(params.out, sep="\t",
		header=True, quoting=False,
		index=False)




#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hp:c:o:', \
			["help", "out=", "paf=", "coords="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.coords = None
		self.paf = None
		self.out = "out.txt"

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "h" or opt == "help":
				continue
			elif opt == "c" or opt=="coords":
				self.coords=arg
			elif opt == "p" or opt=="paf":
				self.paf=arg
			elif opt =="o" or opt=="out":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.paf:
			self.display_help("No paf provided.")
		if not self.coords:
			self.display_help("No coordinates provided.")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nliftoverFromPafscaff.py\n")
		print("Author: Tyler K Chafin, University of Colorado")
		print ("Contact: tyler.chafin@colorado.edu")
		print ("Description: Converts a given set of coordinates (e.g., from a VCF file) to a new coordinate system, as mapped by pafscaff")
		print("""
Arguments:
	-h, --help	: Display help menu
	-p,--paf	: Path to pafscaff fasta file (can be just headers)
	-c,--coords	: Tab-delimited table in the format: scaffold_name "\t" position
	-o, --out	: Output file name [default=out.tsv]
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
