#!/usr/bin/python

import re
import sys
import os
import getopt
from itertools import product
import Bio
from Bio import AlignIO

def main():
	params = parseArgs()

	if params.infile:
		file_object = open("out.fasta", "w")
		d=dict()
		for seq in read_fasta(params.infile):
			loc_name = getLocName(seq[0])
			if loc_name not in d:
				d[loc_name] = list()
			for i in expandAmbiquousDNA(seq[1]):
				d[loc_name].append(i)
		for key, value in d.items():
			aln = Bio.Align.MultipleSeqAlignment([])
			if len(value)>1:
				for seq in value:
					aln.add_sequence("NA",seq)
				locus = consensAlign(aln, threshold=0.1, mask=0.1)
				tr_counts = seqCounterSimple(simplifySeq(locus.conSequence))
				if tr_counts["*"] < params.maxS:
					header = ">" + key + "\n"
					data = locus.conSequence + "\n"
					file_object.write(header)
					file_object.write(data)
			else:
				header = ">" + key + "\n"
				data = seq + "\n"
				file_object.write(header)
				file_object.write(data)

		file_object.close()

	else:
		sys.exit("No input provided.")


#Function to simplify a sequence
def simplifySeq(seq):
	temp = re.sub('[ACGT]', '', (seq).upper())
	return temp.translate(str.maketrans("RYSWKMBDHV", "**********"))

#Returns dict of character counts from a simplified consensus sequence
def seqCounterSimple(seq):
	d = {}
	d = {
		'N':0,
		'-':0,
		'*':0
	}
	for c in seq:
		if c in d:
			d[c] += 1
	return d

def getLocName(header):
	l = header.split("_")
	return(str(l[1]) + "_" + str(l[2]))


class consensAlign():
	'Consensus alignment object'
	#Default constructor
	def __init__(self, alignment, threshold, mask):
		self.alnVars = []
		self.conSequence = make_consensus(alignment, threshold, mask)
		self.alnVars = get_vars(self.conSequence)

class variablePosition():
	'Object to hold information about a variable position'
	#Default constructor
	def __init__(self, pos=None, val=None):
		self.position = pos
		self.value = val.upper()


	@classmethod
	def from_list(cls, data):
		pos = int(data[0])
		val = str(data[1])
		new = cls(pos, val)
		return new



######################## STATIC FUNCTIONS ##############################
#Function to count number of lower case in a string
def n_lower_chars(string):
    return sum(1 for c in string if c.islower())
#Function to return sorted unique string from list of chars
def listToSortUniqueString(l):
	sl = sorted(set(l))
	return(str(''.join(sl)))

#Less shitty consensus function than BioPython has..
#From an AlignIO alignment object
def make_consensus(alignment, threshold=0.1, mask=0.1):
	aln_depth = len(alignment)
	#If only one sequence in alignment, return that seq as consensus
	if aln_depth == 1:
		return alignment[:,0]
	aln_len = alignment.get_alignment_length()
	consensus="" #consensus string to build
	#For each column
	for i in range(aln_len):
		#For each character in column, grab occurences
		col_len = len((alignment[:,i]))
		nuc_count = dict() #dictionary to track occurences
		nuc_types = 0 #track number of things we found
		ismask = 0 #Track if we should mask this column
		#Check number of masked bases
		nlower = n_lower_chars(alignment[:,i])
		prop_mask = float(nlower/aln_depth)
		if prop_mask > mask:
			ismask = 1

		for c in ((alignment[:,i]).upper()):
			for nuc in get_iupac(c):
				nuc_count[nuc] = nuc_count.get(nuc,0)+1
				#print(nuc, end='', flush=True)
		nucs= []
		add = 0
		for key_raw in nuc_count:
			key = key_raw
			if ismask: #If masked above threshold, set to lower case
				key = key_raw.lower()
			#print(key, " is ", nuc_count[key])
			#If only one nuc type, keep it
			if nuc_types == 1:
				consensus+=str(key)
				add=1
				break
			#If N or gap, call consensus N or gap if above threshold
			elif key in ("N", "n", "-"):
				temp = key
				if key == "n":
					temp = "N"
				if float((nuc_count[temp])/aln_depth) >= threshold:
					#print("Found ", nuc_count[key], key, "'s in alignment")
					#print((nuc_count[key])/aln_depth)
					consensus+=str(key)
					add=1
					break
			else:
				nucs.append(str(key))
		if add == 0:
			temp = listToSortUniqueString(nucs)
			consensus+=reverse_iupac_case(temp)
	return(consensus)

#Function to get a list of variablePositions
def get_vars(con):
	#print("Parsing: ", con)
	var_objects = [] #empty list for holding var objects
	#For each position
	for i in range(len(con)):
		#print(i, " is ", con[i])
		#If not monomorphic
		if con[i].upper() not in {"A", "G", "T", "C", "X"}:
			var_objects.append(variablePosition(i, con[i].upper()))
			#continue
		#else:
		'''
		#For each sequence
		for c in range(len(aln[:,i])):
			#print(aln[c,i], end='', flush=True)
			ref = con[i].upper()
			var = aln[c,i].upper()
			#ref.upper()
			#var.upper()
			#print("Var is ",var, " and Ref is ", ref)
			if var == ref:
			#if var == "-" and ref == "-":
				continue
			#elif var == "N" and ref == "N":
				continue
			else:
				#print(var, end='', flush=True)
				#print(aln[c].id, " has ", aln[c,i], " at pos ", i)
				var_objects.append(variablePosition(aln[c].id, i, aln[c,i]))
		'''
	return var_objects

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac(char):
	iupac = {
		"A"	: ["A"],
		"G"	: ["G"],
		"C"	: ["C"],
		"T"	: ["T"],
		"N"	: ["N"],
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
	return iupac[char]

#Function to translate a string of bases to an iupac ambiguity code
def reverse_iupac(char):
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

#Function to translate a string of bases to an iupac ambiguity code, retains case
def reverse_iupac_case(char):
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
		'ACGT':'N',
		'a':'a',
		'n':'n',
		'c':'c',
		'g':'g',
		't':'t',
		'ag':'r',
		'ct':'y',
		'ac':'m',
		'gt':'k',
		'at':'w',
		'cg':'s',
		'cgt':'b',
		'agt':'d',
		'act':'h',
		'acg':'v',
		'acgt':'n'
	}
	return iupac[char]



#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'i:m:h', \
			["in=","max="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.infile=None
		self.maxS=1

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
			if opt in ('i', "in"):
				self.infile = arg
			elif opt in ('m', "max"):
				self.maxS = int(arg)
			elif opt in ('h', 'help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.infile:
			sys.exit("Error: Input file required")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ncollapse_baits.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-i <baitsfile> \n")
		print ("Description: Collapses and filters (by SNP count) baits sampled by BaitsTools pyrad2baits")

		print("""
	Input options:
		-i,--in	: Input .fa file. Assumes loci consist of sequential lines
		-m,--max	: Maximum allowed SNPs per consensus bait [default=1]
		-h,--help	: Displays help menu""")
		print()
		sys.exit()

#Function to split character to IUPAC codes, assuing diploidy
def get_iupac_caseless(char):
	lower = False
	if char.islower():
		lower = True
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
	if lower:
		ret = [c.lower() for c in ret]
	return ret

#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
def read_fasta(fas):
	if not fileCheck(fas):
		raise FileNotFoundError("Fatal exception, file %s not found."%fas)

	fh = open(fas)
	try:
		with fh as file_object:
			contig = ""
			seq = ""
			for line in file_object:
				line = line.strip()
				if not line:
					continue
				line = line.replace(" ","")
				#print(line)
				if line[0] == ">": #Found a header line
					#If we already loaded a contig, yield that contig and
					#start loading a new one
					if contig:
						yield([contig,seq]) #yield
						contig = "" #reset contig and seq
						seq = ""
					contig = (line.replace(">",""))
				else:
					seq += line
		#Iyield last sequence, if it has both a header and sequence
		if contig and seq:
			yield([contig,seq])
	finally:
		fh.close()

#Function to check if a file path is valid
def fileCheck(f):
	return (os.path.isfile(f))

#Function to expand ambiguous sequences
#Generator function
def expandAmbiquousDNA(sequence):
	for i in product(*[get_iupac_caseless(j) for j in sequence]):
		yield("".join(i))


#Call main function
if __name__ == '__main__':
    main()
