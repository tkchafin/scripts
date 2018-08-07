#!/usr/bin/python

import sys
import os
import getopt
import collections

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
		print("ERROR: Popmap file must be provided.")
		sys.exit(1)

	print("Writing new FASTA files...")
	#For each pop, write a new FASTA
	seen = list(seqs.keys())
	pops = make2Dpopmap(pop_assign)
	for pop in pops.keys():
		fas = str(pop) + ".fasta"
		with open(fas, 'w') as fh:
			try:
				print(fas + "....")
				for sample in pops[pop]:
					if sample in seen:
						to_write = ">" + str(sample) + "\n" + seqs[sample] + "\n"
						fh.write(to_write)
					else:
						print("Sample not found in FASTA:",sample)
			except IOError as e:
				print("Could not read file:",e)
				sys.exit(1)
			except Exception as e:
				print("Unexpected error:",e)
				sys.exit(1)
			finally:
				fh.close()

#Makes a dict of lists from a popmap
def make2Dpopmap(p):
	ret = dict()
	for s in p:
		if p[s] not in ret:
			ret[p[s]] = list()
		ret[p[s]].append(s)
	return(ret)



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
			options, remainder = getopt.getopt(sys.argv[1:], 'f:p:h', \
			["ppmap=","fasta=","help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.popmap=None
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
			if opt in ('p', 'popmap'):
				self.popmap = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('f', 'fasta'):
				self.fasta = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.popmap:
			self.display_help("Error: Need popmap")
		if not self.fasta:
			self.display_help("Error: Need fasta")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nsplitFastaPops.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-f <fasta> -p <popmap> \n")
		print ("Description: Splits a FASTA file into 1 file per population (pops from tab-delited popmap)")

		print("""
	Arguments:
		-i,--input	: FASTA file
		-p,--popmap	: Tab-delimited population map (Sample \\t PopID)
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
