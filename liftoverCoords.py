#!/usr/bin/python

import sys
import os
import getopt
import pyliftover
import csv
import pandas as pd


def main():
	params = parseArgs()
	if params.liftover:
		lo = pyliftover.LiftOver(params.liftover)
		if params.table:
			tab=pd.read_csv(params.table, sep="\t")
			print("Read table:")
			print(tab)
			def convert(row):
				name="chr"+row[params.chrom]
				ret=lo.convert_coordinate(name, row[params.bp])
				return(int(ret[0][1]))
				
			tab[params.ocol] = tab.apply(convert,axis = 1)
			print("Writing the output table:")
			print(tab)
			tab.to_csv(params.oname, sep="\t", index=False)
			
			if params.marey:
				marey=make_marey(tab, params.chrom, params.ocol)
				print("Created the following Marey Map input:")
				print(marey)
				mout=params.oname+"_mmap.txt"
				marey.to_csv(mout, sep=" ", quoting=csv.QUOTE_NONNUMERIC, index=False)
			
		else:
			params.display_help("Error: No table provided")
	else:
		params.display_help("Error: No liftover file provided")

#function writes a spoof marey map file from a table of :
#chr \t bp \t cM \t liftover.bp
def make_marey(table, chrom, bp):
	ret=pd.DataFrame()
	ret["map"] = "chr"+table[chrom].astype(str)
	ret["set"] = "fakeset"
	ret["mkr"] = "fakemarker"
	ret["phys"] = table[bp].astype(int)
	ret["gen"] = table["cM"].astype(float)
	return(ret)
		

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hf:t:p:c:n:o:m', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.table = None
		self.liftover = None
		self.chrom = "chr"
		self.bp = "bp"
		self.ocol = "liftover.bp"
		self.oname = None
		self.marey=False

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
			elif opt == "f":
				self.liftover=arg
			elif opt == "t":
				self.table = arg
			elif opt == "p":
				self.bp=arg
			elif opt == "c":
				self.chrom=str(arg)
			elif opt == "n":
				self.ocol=str(arg)
			elif opt == "o":
				self.oname=str(arg)
			elif opt == "m":
				self.marey=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.liftover or not self.table:
			self.display_help("No files provided.")
		self.oname=self.table + ".liftover"


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nliftoverCoords.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Converts a table of physical positions from one genome assembly to another given an \".over.chain.gz\" database")
		print("""
Arguments: 
	-h, --help	: Display help menu
	-f			: Path to .over.chain.gz file
	-t			: Tab-delimited table including coordinates
	-p			: Column name in table containing the physical (bp) coordinates
				   [default = \"bp\"]
	-c			: Column name in table containing the chromosome names
				   [default = \"chr\"]
	-n			: Output column name for new table
				   [default = \"liftover.bp\"]
	-o			: Output file name
				   [default = \"<infile>.liftover\"]
	-m			: (Boolean) Additionally output Marey-Map input file

	NOTE: Chromosomes should be named e.g. as \"chr1\" or \"chrX\" in the
			      .over.chain.gz file, but without the \"chr\" in the table file """)
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
