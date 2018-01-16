#!/usr/bin/python

import re
import sys
import os
import getopt
from itertools import product

def main():
	params = parseArgs()

	if params.fasta:
		file_object = open("out.fasta", "w")
		for seq in read_fasta(params.fasta):
			count = 1
			for i in expandAmbiquousDNA(seq[1]):
				header = ">" + str(seq[0]) + "." + str(count)+ "\n"
				sequence = str(i) + "\n"
				file_object.write(header)
				file_object.write(sequence)
				count += 1
		file_object.close()
	elif params.seq:
		for i in expandAmbiquousDNA(params.seq):
			print(i)
	else:
		sys.exit("No input provided.")


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 's:f:h', \
			["seq=","fasta=","help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.seq=None
		self.fasta=None

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
			if opt in ('s', 'seq'):
				self.seq = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('f','fasta'):
				self.fasta = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if self.seq and self.fasta:
			sys.exit("Error: Input either -s, or -f. Not both.")
		if not self.seq and not self.fasta:
			sys.exit("Error: Input either -s, or -f.")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nexpandSeq.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-s AGTGATAGTAGTGRRTGAYAGAGT \n")
		print ("Description: expandSeq.py expands DNA sequences with ambiguities to a list of all possible variants.")

		print("""
	Input options:
		-s,--seq	: Sequence string to expand (results output to stdout)
			   or
		-f,--fasta	: You can also specify a FASTA file. Results will be output as FASTA.
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
