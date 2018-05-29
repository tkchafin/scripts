#!/usr/bin/python

import re
import sys
import os
import getopt
import operator
import random

def main():
	params = parseArgs()

	if params.phylip:
		#Get sequences as dict of lists
		seqs = readPhylip(params.phylip)

		#get list of columns and list of samplenames
		alen = getSeqLen(seqs)
		columns = [[]for i in range(alen)]
		names = list()
		for key, value in seqs.items():
			names.append(key)
			for i, nuc in enumerate(value):
				columns[i].append(nuc)

		#For each column, delete those which are not bi-allelic
		dels=list()
		for i, col in enumerate(columns):
			if not isBiallelic(col):
				dels.append(i)
				#print(i,"not biallelic:",col)

		print("Deleting",len(dels),"non-biallelic columns.")
		for col in sorted(dels,reverse=True): #reverse sorted so subsequent deletes aren't thrown off
			#print(col,":",columns[col])
			del columns[col]

		#Then, convert to 012 format
		print("Converting to 012 format...")
		formatted = [[]for i in range(alen-len(dels))]

		for i, col in enumerate(columns):
			#print(col)
			#print(nucs2numeric(col))
			if params.nohet:
				formatted[i] = nucs2numericNohet(col)
			else:
				formatted[i] = nucs2numeric(col)
			#sys.exit()

		final_data = dict()
		for i, samp in enumerate(names):
			seqs = list()
			for k,nuc in enumerate(formatted):
				seqs.append(nuc[i])
			final_data[samp] = "".join(seqs)

		print("Writing NEXUS output file...")
		dict2nexus(params.out, final_data)


	else:
		print("No input provided.")
		sys.exit(1)

#Function takes biallelic list of nucleotides and converts to numeric
#0 = major allele
#1 = minor allele
#2 = het
#? = - or N
def nucs2numeric(nucs):
	if isBiallelic(nucs):
		#print(nucs)
		ret = list()
		counts = {"A":0, "G":0, "C":0, "T":0}
		#find major allele
		for nuc in nucs:
			if nuc not in ("-", "N"):
				for exp in get_iupac_caseless(nuc):
					counts[exp] += 1
		#sort dict, to list of tuples (b/c dicts are orderless, can't keep as dict)
		sorted_x = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
		majA = sorted_x[0][0]
		minA = sorted_x[1][0]
		het = reverse_iupac(''.join(sorted(set([majA, minA])))) #get het code
		#print(majA, minA, het)

		for nuc in nucs:
			nuc = nuc.upper()
			if nuc == majA:
				ret.append("0")
			elif nuc == minA:
				ret.append("1")
			elif nuc == het:
				ret.append("2")
			elif nuc == "-":
				ret.append("-")
			else:
				ret.append("?")

		return(ret)
	else:
		print("Warning: Data is not biallelic:",nucs)
		return(None)

#Function takes biallelic list of nucleotides and converts to numeric
#0 = major allele
#1 = minor allele
#2: Randomly samples heterozygous sites as 0 or 1
def nucs2numericNohet(nucs):
	if isBiallelic(nucs):
		#print(nucs)
		ret = list()
		counts = {"A":0, "G":0, "C":0, "T":0}
		#find major allele
		for nuc in nucs:
			if nuc not in ("-", "N"):
				for exp in get_iupac_caseless(nuc):
					counts[exp] += 1
		#sort dict, to list of tuples (b/c dicts are orderless, can't keep as dict)
		sorted_x = sorted(counts.items(), key=operator.itemgetter(1), reverse=True)
		majA = sorted_x[0][0]
		minA = sorted_x[1][0]
		het = reverse_iupac(''.join(sorted(set([majA, minA])))) #get het code
		#print(majA, minA, het)

		for nuc in nucs:
			nuc = nuc.upper()
			if nuc == majA:
				ret.append("0")
			elif nuc == minA:
				ret.append("1")
			elif nuc == het:
				ret.append(random.randint(0,1))
			elif nuc == "-":
				ret.append("-")
			else:
				ret.append("?")

		return(ret)
	else:
		print("Warning: Data is not biallelic:",nucs)
		return(None)

#Function to translate a string of bases to an iupac ambiguity code
def reverse_iupac(char):
	char = char.upper()
	iupac = {
		'A':'A',
		'N':'N',
		'-':'-',
		'C':'C',
		'G':'G',
		'T':'T',
		'AG':'R',
		'CT':'Y',
		'AC':'M',
		'GT':'K',
		'AT':'W',
		'CG':'S',
		'CGT':'B',
		'AGT':'D',
		'ACT':'H',
		'ACG':'V',
		'ACGT':'N'
	}
	return iupac[char]

#Function takes a list of nucleotides, and returns True if the column is biallelic
#ignores gaps and Ns
#expands uipac codes using a call to external function
def isBiallelic(nucs):
	expanded = list()
	for nuc in nucs:
		if nuc not in ("-", "N"):
			for exp in get_iupac_caseless(nuc):
				expanded.append(exp)
	uniq_sort = sorted(set(expanded))
	if len(uniq_sort) != 2:
		#print(nucs)
		#print(uniq_sort, len(uniq_sort))
		return(False)
	else:
		return(True)

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
	if char.islower():
		char = char.upper()
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["A", "C", "G", "T"],
		"-"	: ["-"],
		"R"	: ["A","G"],
		"Y"	: ["C","T"],
		"S"	: ["G","C"],
		"W"	: ["A","T"],
		"K"	: ["G","T"],
		"M"	: ["A","C"],
		"B"	: ["C","G","T"],
		"D"	: ["A","G","T"],
		"H"	: ["A","C","T"],
		"V"	: ["A","C","G"]
	}
	ret = iupac[char]
	return ret

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


#Function to write an alignment as DICT to NEXUS
def dict2nexus(nex, aln):
	with open(nex, 'w') as fh:
		try:
			slen = getSeqLen(aln)
			header = "#NEXUS\n\nBegin data;\nDimensions ntax=" + str(len(aln)) + " nchar=" + str(slen) + ";\n"
			header = header + "Format datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n\n"
			fh.write(header)
			for seq in aln:
				sline = str(seq) + " " + aln[seq] + "\n"
				fh.write(sline)
			last = ";\nEnd;\n"
			fh.write(last)
		except IOError:
			print("Could not read file ",nex)
			sys.exit(1)
		finally:
			fh.close()

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


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'p:ho:n', \
			["phylip=","phy=","out=","nohet"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.phylip=None
		self.out="out.nex"
		self.nohet=False

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
			if opt in ('p', 'phylip', 'phy'):
				self.phylip = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ('n','nohet'):
				self.nohet=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.phylip:
			self.display_help("Error: Missing required phylip file (-p, --phylip)")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nphylip2biNumNex.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-p /path/to/phylip \n")
		print ("Description: Converts PHYLIP file to NEXUS file of only bi-allelic markers, coded with 012. As inputs for PhyloNetworks MLE_biMarkers or SNAPP")

		print("""
	Arguments:
		-p,--popmap	: Path to tab-delimited population map
		-o,--out	: Output file name <default = out.nex>
		-n,--nohet	: Randomly sample one allele from all heterozygous sites
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
