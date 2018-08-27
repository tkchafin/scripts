#!/usr/bin/python

import re
import sys
import os
import getopt
import operator
import collections
import copy

def main():
	params = parseArgs()

	if params.phylip:
		#Get sequences as dict of lists
		seqs = readPhylip(params.phylip)
	else:
		print("No input provided.")
		sys.exit(1)

	pop_assign = dict()
	#parse popmap file for dictionary of sample assignments
	if params.popmap:
		print("Parsing popmap file...")
		pop_assign = parsePopmap(params.popmap)
	else:
		print("ERROR: Popmap file must be provided.")
		sys.exit(1)

	if seqs and pop_assign:

		#Remove samples from pop_assign that do not have data
		pop_assign = cleanPopmap(pop_assign, seqs.keys())

		#Make dict of dicts that splits by population, only retaining pops/samples from the popmap.
		#Get unique pop names using one of the worst python lines ever written
		pops = dict()
		for k in set(pop_assign.values()):
			pops[k] = dict()

		#Remove pops listed as excluded
		if params.exclude:
			print("Excluding populations:", ", ".join(params.exclude))
			for exc in params.exclude:
				if exc in pops:
					del pops[exc]
		if params.include:
			print("Only keeping populations:", ", ".join(params.include))
			for pop in list(pops):
				if pop not in params.include:
					del pops[pop]

		#make sure we didn't throw out all populations...
		if len(list(pops)) < 1:
			print("Oops! No populations remaining. Check that popmap sample names match those in your data file, or that selections using --include or --exclude are correct! :)")
			sys.exit(1)

		alen = getSeqLen(seqs)
		inum = 0
		for assigned in pop_assign:
			if pop_assign[assigned] in pops:
				pops[pop_assign[assigned]][assigned] = seqs[assigned]
				inum+=1
		seqs.clear()

		#Make 2D list to remove columns failing the globalN filter
		bad_columns = list() #list of column numbers to delete

		#For each pop dict, make 2D list to remove columns failing popN filter
		print("Found",alen,"nucleotide columns in the dataset!")
		columns = [[]for i in range(alen)] #2D array of global data
		for pop, data in pops.items():
			for sample, sequence in data.items():
				for i, nuc in enumerate(sequence):
					columns[i].append(nuc)

		#Write new ordered output and phylip
		print("Writing outputs...")
		phy = params.out + ".phy"
		omap = params.out + ".map"

		pfh = open(phy, "w")
		mfh = open(omap, "w")

		header = str(inum) + "\t" + str(alen) + "\n"
		pfh.write(header)

		for pop in sorted(pops):
			for ind, data in pops[pop].items():
				indline = str(ind) + "\t" + "".join(data) + "\n"
				pfh.write(indline)

				mapline = str(pop) + "\n"
				mfh.write(mapline)
		pfh.close()
		mfh.close()

		print("Done!\n")

#Goes through a dict of sequences and get the alignment length
def getSeqLen(aln):
	length = None
	for key in aln:
		if not length:
			length = len(aln[key])
		else:
			if length != len(aln[key]):
				print("getSeqLen: Alignment contains sequences of multiple lengths.")
	return(length)

#function reads a tab-delimited popmap file and return dictionary of assignments
def parsePopmap(popmap):
	ret = dict()
	if os.path.exists(popmap):
		with open(popmap, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					else:
						stuff = line.split()
						ret[stuff[0]] = stuff[1]
				return(ret)
			except IOError:
				print("Could not read file ",pairs)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%popmap)

#Function to remove samples from a popmap dict, given a list of valid samples (e.g. those to retain)
def cleanPopmap(popmap, names):
	ret = copy.deepcopy(popmap)
	to_remove = list()
	for ind in popmap:
		if ind not in names:
			to_remove.append(ind)
	for rem in sorted(to_remove, reverse=True):
		del ret[rem]

	return(ret)

#Function to read a phylip file. Returns dict (key=sample) of lists (sequences divided by site)
def readPhylip(phy):
	if os.path.exists(phy):
		with open(phy, 'r') as fh:
			try:
				num=0
				ret = dict()
				for line in fh:
					line = line.strip()
					if not line:
						continue
					num += 1
					if num == 1:
						continue
					arr = line.split()
					ret[arr[0]] = list(arr[1])
				return(ret)
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'p:i:ho:X:I:', \
			["input=","phylip=","phy=","out=","popmap=","maxN=",
			"popN=","exclude=","include="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.phylip=None
		self.popmap=None
		self.out="out"
		self.exclude = list()
		self.include = list()


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
			if opt in ('i', 'phylip', 'input','phy'):
				self.phylip = arg
			elif opt in ('p', 'popmap'):
				self.popmap = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ('X', 'exclude'):
				self.exclude = arg.split(",")
			elif opt in ('I','include'):
				self.include = arg.split(",")
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.phylip and not self.fasta:
			self.display_help("Error: Missing required alignment file (--fasta or --input)")
		if not self.popmap:
			self.display_help("Error: Missing required popmap file (-p, --popmap)")
		if self.include and self.exclude:
			self.display_help("Don't use both --include and --exclude.")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nmakeHyde.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-i /path/to/phylip -i /path/to/popmap\n")
		print ("Description: Making inputs for HyDe and filtering populations for inclusion/exclusion")

		print("""
	Arguments:
		INPUT FILES [REQUIRED]
		-i,--input	: Input file as PHYLIP
		-p,--popmap	: Tab-delimited population map

		PARAMETERS [OPTIONAL]
		-o,--out	: Output file name <default = out.nex>
		-X,--exclude: List of pops to exclude (format: -x "Pop1,Pop2,Sample4...")
		-I,--include: List of pops to include (removing all others)
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
