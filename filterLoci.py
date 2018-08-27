#!/usr/bin/python

import re
import sys
import os
import getopt


def main():
	params = parseArgs()

	if params.input:

		allLoci = list()
		locNum = 1
		with open(params.input, "r") as fh:
			try:
				aln_d = dict()
				numSamp = int(0)
				numPIS = int(0)

				passed=0
				sampFail=0
				varsFail=0

				for line in fh:
					line = line.strip()
					if not line:
						continue
					if line[0] == ">":
						allLoci.append(line)
						line = line.replace(">","")
						stuff = line.split()
						identifier = stuff[0].replace(" ", "")
						seqs = stuff[1].replace(" ", "")
						aln_d[identifier] = seqs
						numSamp += 1
					elif line[0] == "/": #alignment end. grab numPIS and parse locus
						numPIS = line.count("*")
						allLoci.append(line)
						locNum += 1
						if numPIS >= params.pis and numSamp >= params.samples:
							passed += 1
							#Write to desired location
							if params.fas:
								f = str(params.out) + "_" + str(locNum) + ".fasta"
								writeFasta(aln_d, f)

							if params.nex:
								n = str(params.out) + "_" + str(locNum) + ".nex"
								dict2nexus(n, aln_d)
							numSamp=0
							numPIS=0
						else:
							if numPIS < params.pis:
								varsFail += 1
							if numSamp < params.samples:
								sampFail += 1

			except IOError as e:
				print("Couldn't read file %s: %s" %(params.input,e))
				sys.exit(1)
			except Exception as e:
				print("Unexpected error reading file %s: %s" %(params.input, e))
				sys.exit(1)
			finally:
				fh.close()

		if params.loci:
			l = str(params.out) + "_filtered.loci"
			writeStuff(allLoci, l)
		print("Number of loci passing filters: %s (of %s total)" %(passed, locNum))
		print("\t",sampFail,"loci failed for too few individuals.")
		print("\t",varsFail,"loci failed for too few parsimony-informative sites\n")

	else:
		print("No input provided.")
		sys.exit(1)

#function writes everything in list to a file
def writeStuff(stuff, f):
	with open(f, "w") as fh:
		try:
			for line in stuff:
				out = str(line) + "\n"
				fh.write(out)
		except IOError as e:
			print("Couldn't write to file %s: %s" %(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error writing to file %s: %s" %(f, e))
			sys.exit(1)
		finally:
			fh.close()

#function writes sequences in dictionary d to fasta-file f
def writeFasta(d, f):
	with open(f, "w") as ffh:
		try:
			for samp, seq in d.items():
				l = ">" + str(samp) + "\n" + seq + "\n"
				ffh.write(l)
		except IOError as e:
			print("Couldn't write to file %s: %s" %(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error writing to file %s: %s" %(f, e))
			sys.exit(1)
		finally:
			ffh.close()

#Function to write an alignment as DICT to NEXUS
def dict2nexus(nex, aln):
	with open(nex, 'w') as fh:
		try:
			slen = getSeqLen_min(aln)
			header = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=" + str(len(aln)) + " NCHAR=" + str(slen) + ";\n"
			header = header + "FORMAT DATATYPE=DNA MISSING=N GAP=-;\n\nMATRIX\n"
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
#returns minimum length
def getSeqLen_min(aln):
	length = None
	for key, val in aln.items():
		if not length:
			length = len(val)
		else:
			if length >= len(val):
				length = len(val)
	return(length)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'i:s:p:ho:nfl', \
			["input=", "samples=", "pis=","help","out=", "nex", "fas", "loci"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.input=None
		self.samples=1
		self.pis=1
		self.out="loc"

		self.nex=False
		self.fas=False
		self.loci=False


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
			if opt in ('i', 'input'):
				self.input = arg
			elif opt in ('p', 'pis'):
				self.pis=int(arg)
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ('s', 'samples'):
				self.samples = int(arg)
			elif opt in ('n', 'nex'):
				self.nex=True
			elif opt in ('f', 'fas'):
				self.fas=True
			elif opt in ('l', 'loci'):
				self.loci=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.input:
			self.display_help(".loci file must be provided (-l, --loci)")

		if not self.loci and not self.nex and not self.fas:
			print("Warning: No output type chosen. Are you sure you wanted this?? filterLoci.py will only report number of loci passing filters.")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nfilterLoci.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-l /path/to/loci\n")
		print ("Description: Filters pyRAD .loci file for N SNPs and N samples")

		print("\nProbably will add more functions later...")

		print("""
	Arguments:
		INPUT FILES [REQUIRED]
		-i,--input	: Input file as PHYLIP

		PARAMETERS [OPTIONAL]
		-s,--samples	: Minimum number of samples to retain locus [def=1]
		-p,--pis	: Minimum number of parsimony-informative sites [def=1]
		-n,--nex	: Boolean. Write loci as individual NEXUS files [false]
		-f,--fas	: Boolean. Write loci as individual FASTA files [false]
		-l,--loci	: Boolean. Write loci as new .loci file [false]
		-o,--out	: Output file prefix

""")
		sys.exit()


#Call main function
if __name__ == '__main__':
    main()
