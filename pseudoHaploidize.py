#!/usr/bin/python

import re
import sys
import os
import getopt
import random

def main():
	params = parseArgs()

	seqs = dict() #key=FASTA header; val=sequence

	#Now, get the alignment from the FASTA file
	#note that this works fine with interleaved FASTA
	if params.fasta:
		print('Reading alignment from FASTA...')
		for f in read_fasta(params.fasta):
			seqs[f[0]] = list(f[1])


	#get indices of all multi-allele sites, then randomly resolve each
	mults = ["R", "Y", "S", "W", "K", "M", "D", "H", "B", "V"]

	for key in (seqs.keys()):
		#get indices of multi-allelic sites
		idxs = [i for i, c in enumerate(seqs[key]) if c.upper() in mults]

		#loop through amiguities, replace each with a new one
		for i in idxs:
			#print(seqs[key][i], end=" - ")
			seqs[key][i] = sampleAllele(seqs[key][i])
			#print(seqs[key][i])

	#write new FASTA outputs
	for samp in seqs.keys():
		seqs[samp] = "".join(seqs[samp])
	if (params.split):
		for samp in seqs.keys():
			fname = samp + "_" + params.out
			sd = dict()
			sd[samp] = seqs[samp]
			write_fasta(fname, sd)
	else:
		write_fasta(params.out, seqs)

#Function to write fasta-formatted sequences
def write_fasta(f, aln):

	with open(f, 'w') as fh:
		try:
			for samp in aln.keys():
				ol = ">" + str(samp) + "\n" + str(aln[samp]) + "\n"
				fh.write(ol)
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()

#function to randomly sample an allele given an ambiguity code
def sampleAllele(ch):
	return(random.choice(get_iupac(ch.upper())))

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

#function returns all indices
def find(str, opts):
	return [i for i, ltr in enumerate(s) if ltr == ch]

#Read samples as FASTA. Generator function
def read_fasta(fas):

	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#print(line)
					if line[0] == ">": #Found a header line
						#If we already loaded a contig, yield that contig and
						#start loading a new one
						if contig:
							yield([contig,seq]) #yield
							contig = "" #reset contig and seq
							seq = ""
						split_line = line.split()
						contig = (split_line[0].replace(">",""))
					else:
						seq += line
				#Iyield last sequence, if it has both a header and sequence
				if contig and seq:
					yield([contig,seq])
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
			options, remainder = getopt.getopt(sys.argv[1:], 'f:so:h', \
			["out=" "help", "fasta=", "split"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.fasta=None
		self.out=None
		self.split=False

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
			if opt =="f" or opt=="fasta":
				self.fasta = arg
			elif opt =="o" or opt=="out":
				self.out = arg
			elif opt == "s" or opt == "split":
				self.split=True
			elif opt =="h" or opt == "help":
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasta:
			self.display_help("Must provide FASTA file <-f,--fasta>")

		#get output prefix if not set by user
		if not self.out:
			self.out = os.path.splitext(self.fasta)[0] + '_hap.fasta'

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\npseudoHaploidize.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-f <input.fasta> [-s] [-f example_hap]\n")
		print ("Description: Creates a pseudo-haploid sequence from input fasta, randomly resolving heterozygous sites")

		print("""
	Arguments:
		-f,--fasta	: Input fasta sequence
		-s,--split	: [Boolean] Write outputs each to their own output file
		-o,--out	: Output file name [default=input_hap.fasta or samp_input_hap.fasta if -s]
		-h,--help	: Displays help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
