#!/usr/bin/python

import re
import sys
import os
import getopt
import random

def main():
	params = parseArgs()

	seqs = dict() #key=FASTA header; val=sequence

	#read sequence in
	if params.fasta:
		print('Reading alignment from FASTA...')
		for f in read_fasta(params.fasta):
			seqs[f[0]] = f[1]

		print("Writing new PHYLIP file",params.out)
		write_phylip(params.out, seqs)
	elif params.phylip:
		print('Reading alignment from PHYLIP...')
		for f in read_phylip(params.phylip):
			seqs[f[0]] = f[1]

		print("Writing new FASTA file",params.out)
		write_fasta(params.out, seqs)



#Print dict to phylip file
def write_phylip(p, aln):
	with open(p, 'w') as fh:
		try:
			header = getPhylipHeader(aln) + "\n"
			fh.write(header)

			for sample in aln.keys():
				line = str(sample) + "\t" + "".join(aln[sample]) + "\n"
				fh.write(line)
		except IOError as e:
			print("Could not read file %s: %s"%(p,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(p,e))
			sys.exit(1)
		finally:
			fh.close()

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

#Returns header for Phylip file from a dictionary of samples w/ data
def getPhylipHeader(d):
	numSamp = 0
	numLoci = None
	for sample in d:
		numSamp = numSamp + 1
		if not numLoci:
			numLoci = len(d[sample])
		else:
			if numLoci != len(d[sample]):
				print("getPhylipHeader: Warning: Sequences of unequal length.")
	header = str(numSamp) + " " + str(numLoci)
	if numLoci == 0 or not numLoci:
		print("getPhylipHeader: Warning: No loci in dictionary.")
	if numSamp == 0:
		print("getPhylipHeader: Warning: No samples in dictionary.")
	return(header)

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

#Read samples as PHYLIP. Generator function
def read_phylip(phy):
	if os.path.exists(phy):
		with open(phy, 'r') as fh:
			try:
				num=0
				for line in fh:
					line = line.strip()
					if not line:
						continue
					num += 1
					if num == 1:
						continue
					arr = line.split()
					yield(arr[0], arr[1])
			except IOError:
				print("Could not read file ",phy)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%phy)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'f:p:h', \
			["help", "fasta=", "phy="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.fasta=None
		self.phylip=None
		self.out=None

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
			elif opt =="p" or opt=="phy":
				self.phylip = arg
			elif opt =="h" or opt == "help":
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasta and not self.phylip:
			self.display_help("Must provide either a FASTA or PHYLIP file.")

		if self.fasta and self.phylip:
			self.display_help("Must provide either a FASTA or PHYLIP file.")

		#get output prefix if not set by user
		if self.fasta:
			self.out = os.path.splitext(self.fasta)[0] + '.phylip'
		elif self.phylip:
			self.out = os.path.splitext(self.phylip)[0] + '.fasta'

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfasta2phylip.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "[-f <.fasta>] [-p <.phy>]\n")
		print ("Description: Simple script to convert between FASTA and PHYLIP formats")

		print("""
	Arguments:
		-f,--fasta	: Input FASTA to be converted to PHYLIP
		-p,--phy	: Input PHYLIP to be converted to FASTA
		-h,--help	: Displays help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
