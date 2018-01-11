#!/usr/bin/python

import getopt
import sys
import os

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'x:p:l:o:h', \
			["xml=","phy=","list=","out=","help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.xml=None
		self.phy=None
		self.tax=None
		self.out="out.phy"

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
			if opt in ('x', 'xml'):
				self.xml = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('p','phy'):
				self.phy = arg
			elif opt in ('l','list'):
				self.tax = arg
			elif opt in ('o','out'):
				self.out = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		self.phy or self.display_help("INPUT ERROR: No PHYLIP provided")
		self.tax or self.display_help("INPUT ERROR: No TAXON LIST provided")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nsubsetPhy.py\n")
		print ("Contact:\n\n\tTyler K. Chafin\n\tUniversity of Arkansas\n\ttkchafin@uark.edu\n")
		print ("\nUsage:\n\t", sys.argv[0], "-p </path/to/xml> -l </path/to/.txt\n")
		print ("Description:\n")
		print("\tsubsetPhy.py is a quickly written and shitty script to help manipulate phylip files\n")

		print("""
		Input options:

			-p,--phy	: Phylip file
			-l,--list	: .txt file containing a list of taxa to subset
			-o,--out	: (Optional) output prefix [default:out.xml]
			-h,--help	: Displays help menu""")
		print()
		sys.exit()


################################# MAIN #########################################
params = parseArgs()

#Read TAX LIST into a list
taxlist = list()
fullnames = list()
fh = open(params.tax)
try:
	with fh as file_object:
		for line in file_object:
			line = line.strip()
			if not line:
				continue
			line = line.replace(" ","")
			arr = line.split("_")
			taxlist.append(arr[-1])
			fullnames.append(line)
finally:
	fh.close()

#Read phylip file
data = {}
numSites = None
count = 0
pfh = open(params.phy)
try:
	with pfh as file_object:
		for line in file_object:
			line = line.strip()
			if not line:
				continue
			count += 1
			if count == 1:
				continue
			arr = line.split()
			if arr[0] in taxlist:
				data[fullnames[taxlist.index(arr[0])]] = arr[1]
				if numSites:
					if len(arr[1]) != numSites:
						sys.exit("ERROR: Samples do not have the same sequence length -")
				else:
					numSites = len(arr[1])
finally:
	pfh.close()

#Open output file
ofh = open(params.out, "w")
try:
	with ofh as file_object:
		header = str(len(data)) + " " + str(numSites) + "\n"
		file_object.write(header)
		for key in data:
			out = key + "\t" + data[key] + "\n"
			file_object.write(out)
finally:
	ofh.close()
