#!/usr/bin/python

import sys
import os
import getopt
from itertools import product

def main():
	params = parseArgs()
	
	codon_table = dict()
	
	if params.code is None:
		codon_table = get_standard_code()
	else:
		codon_table = read_code_file(params.code)
	
	#print(codon_table)
	
	nucs = dict()
	current=None
	curr_index=1
	for aa in read_fasta(params.input):
		if current is None:
			current = aa[0]
		elif current != aa[0]:
			current = aa[0]
			curr_index = 1
		for trans in get_all_revtrans(aa[1], codon_table):
			header=str(aa[0]) + "_translation-" + str(curr_index)
			nucs[header] = trans
			curr_index +=1
	
	write_fasta(params.out, nucs)

#generator function
def get_all_revtrans(aa, code):
	possibilities = list()
	for pos in aa:
		possibilities.append(list(code[pos.upper()]))
		
	for nuc in product(*possibilities):
		yield("".join(nuc))

def read_code_file(file):
	d = dict()
	if os.path.exists(file):
		with open(file, 'r') as fh:
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
					if arr[0] not in d:
						d[arr[0].upper()] = list()
					d[arr[0].upper()].append(arr[1].upper())
					
				return(d)
			except IOError:
				print("Could not read file ",file)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%file)
	return(d)


def get_standard_code():
	d = {
		'*' : ['TAA','TAG','TGA'],
		'A' : ['GCA','GCC','GCG','GCT'],
		'C' : ['TGC','TGT'],
		'D' : ['GAC','GAT'],
		'E' : ['GAA','GAG'],
		'F' : ['TTC'],
		'G' : ['GGA','GGC','GGG','GGT'],
		'H' : ['CAC','CAT'],
		'I' : ['ATA','ATC','ATT'],
		'K' : ['AAA','AAG'],
		'L' : ['CTA','CTC','CTG','CTT','TTA','TTG'],
		'M' : ['ATG'],
		'N' : ['AAC','AAT'],
		'P' : ['CCA','CCC','CCG','CCT'],
		'Q' : ['CAA','CAG'],
		'R' : ['AGA','AGG','CGA','CGC','CGG','CGT'],
		'S' : ['AGC','AGT','TCA','TCC','TCG','TCT'],
		'T' : ['ACA','ACC','ACG','ACT'],
		'V' : ['GTA','GTC','GTG','GTT'],
		'W' : ['TGG'],
		'Y' : ['TAC','TAT']
	}
	return(d)

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
			options, remainder = getopt.getopt(sys.argv[1:], 'hi:c:o:', \
			["help", "in=", "code=", "out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.input=None
		self.code=None
		self.out="out.fas"


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
			elif opt=="i" or opt=="in":
				self.input = arg
			elif opt=="c" or opt=="code":
				self.code = arg
			elif opt=="out" or opt=="o":
				self.out = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.input:
			self.display_help("No input provided.")
		
		if not self.code:
			self.display_help("No code provided. Using default.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nrevTransAll.py\n")
		print("Author: Tyler Chafin")
		print ("Contact: tyler.chafin@colorado.edu")
		print ("Description: Gives all possible reverse translations for a amino acid sequence")
		print("""
		-i,--in 	: Input file name (FASTA format)
			format:
			>my_sequence
			MFLIMVVFPTTAASVMMVMMV...
		-c,--code	: Tab-delimited codon table 
			format:
			F	TTT
			F	TTC
			F	TTA
			F	TTG
			L	CTT
			...
			...
			<NOTE: If none supplied, will use 'standard' code>
		-o,--out	: Output file name (default=out.fas)
			format:
			>my_sequence_translation-1
			ATGATGAT...
			>my_sequence_translation-2
			ATGATCAT...
			...
			...
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
	main()
