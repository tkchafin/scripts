#!/usr/bin/python

import re
import sys
import os
import getopt

def main():
	params = parseArgs()

	pop_assign = dict()
	seqs = dict()

	#parse popmap file for dictionary of sample assignments
	if params.popmap:
		print("Parsing popmap file...")
		pop_assign = parsePopmap(params.popmap)
	else:
		print("ERROR: Popmap file must be provided.")
		sys.exit(1)

	#Now, get the alignment from the FASTA file (as another dict)
	if params.fasta:
		print('Reading alignment from FASTA...')
		for f in read_fasta(params.fasta):
			seqs[f[0]] = f[1]
	else:
		if not params.nex:
			print("ERROR: Fasta or Nexus file must be provided.")
			sys.exit(1)

	#Write alignment to NEXUS
	if params.nex:
		#Get sample names, check if they're all in the popmap
		seqnames = getSamplesNexus(params.nex)
		print("Checking that samples in NEXUS match POPMAP...")
		#validatePopmap(seqnames, pop_assign)
	else:
		#Check fasta samples are in popmap
		print("Writing alignment to NEXUS file...")
		dict2nexus(params.out, seqs)

	#Get unique from list, enumerate, convert to dict
	#One of the most horrible lines of Python I've ever written :)
	pops = dict(enumerate(list(set(list(pop_assign.values())))))
	with open(params.out, 'a') as fh:
		try:
			print("Appending PopArt block to NEXUS file...")
			#build header lines
			header = "BEGIN TRAITS;\n\tDIMENSIONS NTRAITS=" + str(len(pops.keys())) + ";\n\tFormat labels=yes missing=? separator=Comma;" + "\n\tTraitLabels"
			for p in pops:
				header = header + " " + str(pops[p])
			header = header + ";\nMATRIX\n"
			fh.write(header)

			#Write poptraits for each sample...
			for ind in pop_assign:
				if ind in seqs.keys():
					iline = ind + " "
					for i in pops:
						if pops[i] == pop_assign[ind]:
							iline = iline + "1"
						else:
							iline = iline + "0"
						if int(i) == len(pops.keys())-1:
							iline = iline + "\n"
						else:
							iline = iline + ","
					fh.write(iline)
			#write end of block
			end = ";\nEND;\n"
			fh.write(end)
		except IOError:
			print("Could not append to file ",params.out)
			sys.exit(1)
		finally:
			fh.close()




#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
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

#Function to read NEXUS and get list of sample names
def getSamplesNexus(nex):
	if os.path.exists(nex):
		with open(nex, 'r') as fh:
			try:
				ret = list()
				found = False
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#print(line)
					if line[0:6].upper() == "MATRIX": #Found a beginning of samples
						found = True
						continue

					if found:
						if line[0] == ";":
							break
						stuff = line.split()
						ret.append(stuff[0])
					else:
						continue
				return(ret)

			except IOError:
				print("Could not read file ",nex)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%nex)

#Function to check that list of sample names and popmap entries match
def validatePopmap(samples, popmap):
	print(samples)
	print(popmap)
	for samp in samples:
		if samp in popmap:
			print("Warning: Sample %s not found in popmap!"%samp)
	for key in popmap:
		if key not in samples:
			print("Warning: Sample %s found in popmap has no data!"%key)

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

#Function to write an alignment as DICT to NEXUS
def dict2nexus(nex, aln):
	with open(nex, 'w') as fh:
		try:
			slen = getSeqLen(aln)
			header = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=" + str(len(aln)) + " NCHAR=" + str(slen) + ";\n"
			header = header + "FORMAT DATATYPE=DNA MISSING=? GAP=-;\n\nMATRIX\n"
			fh.write(header)
			for seq in aln:
				sline = str(seq) + "\t" + aln[seq] + "\n"
				fh.write(sline)
			last = ";\nEND;\n"
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
			options, remainder = getopt.getopt(sys.argv[1:], 'p:f:hn:o:', \
			["popmap=","help","fasta=","nex=","out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.popmap=None
		self.out=None
		self.fasta=None
		self.nex=None

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
			if opt in ('p', 'popmap'):
				self.popmap = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ("f","fasta"):
				self.fasta = arg
			elif opt in ("n","nex"):
				self.nex = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.popmap:
			self.display_help("Error: Missing required popmap file (-p,--popmap).")
		if not self.fasta:
			if not self.nex:
				self.display_help("Error: Either FASTA or NEXUS must be provided.")

		if self.nex:
			self.out = self.nex
		else:
			if self.out:
				self.out = self.out + ".nex"
			else:
				self.out = "out.nex"


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nmakePopMap.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-p </path/to/popmap> -f </path/to/fasta>\n")
		print ("Description: Creates NEXUS input for PopArt, given FASTA aligment and tab-delimited population map.")

		print("""
	Arguments:
		-p,--popmap	: Path to tab-delimited population map
		-o,--out	: Prefix for output file <default = ./out>
		-f,--fasta	: Path to FASTA-formatted haplotype sequences
			--Note for diplotypes:
				Paired haplotypes should be formatted as SampleName_A and _B
		-n,--nex	: Optionally, can provide a NEXUS file to append to.
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
