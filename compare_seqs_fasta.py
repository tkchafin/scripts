#!/usr/bin/python

import sys
import os
import getopt
from textwrap import wrap

def main():
	params = parseArgs()

	seqs=dict()
	seqlen=0
	for s in read_fasta(params.infile):
		seqs[s[0]]=s[1]
		seqlen=len(s[1])

	ann=""
	for i in range(seqlen):
		seen=list()
		for s in seqs:
			seen.append(seqs[s][i].lower())
		vars=set(seen)
		# add conditions here if you want to change output e.g. if there are gaps
		if len(vars)==1:
			ann = ann + "-"
		else:
			ann = ann + "*"
	seqs["annotation"] = ann

	write_fasta(params.out, seqs)

#Function to write fasta-formatted sequences
def write_fasta(f, aln, width=None):
	with open(f, 'w') as fh:
		try:
			for samp in aln.keys():
				if width:
					ol = ">" + str(samp) + "\n"
					chunks=wrap(aln[samp], width=width, break_on_hyphens=False, drop_whitespace=False)
					for chunk in chunks:
						ol=ol + str(chunk) + "\n"
				else:
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
			options, remainder = getopt.getopt(sys.argv[1:], 'hi:o:', \
			["help", "infile=","out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.infile=None
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
				self.infile=arg
			elif opt=="out" or opt=="o":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.infile:
			self.display_help("No files provided.")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("Description: Annotate differences in a fasta alignment")
		print("""
		-i,--in		: Input fasta file
		-o,--out	: Output file name (default=out.fas)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
