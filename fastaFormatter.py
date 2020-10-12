#!/usr/bin/python

import sys
import os
import getopt
from textwrap import wrap

def main():
	params = parseArgs()
	
	if params.many2one:
		seqs=dict()
		for f in read_fasta(params.many2one):
			seqs[f[0]] = f[1]
		write_fasta(params.out, seqs)
	elif params.one2many:
		seqs=dict()
		for f in read_fasta(params.one2many):
			seqs[f[0]] = f[1]
		write_fasta(params.out, seqs, params.width)

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
			options, remainder = getopt.getopt(sys.argv[1:], 'h1:M:w:o:', \
			["help", "one2many=","many2one=","width=","out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.one2many=None
		self.many2one=None
		self.width=60
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
			elif opt=="one2many" or opt=="1":
				self.one2many=arg
			elif opt=="many2one" or opt=="M":
				self.many2one=arg
			elif opt=="width" or opt=="w":
				self.width=int(arg)
			elif opt=="out" or opt=="o":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.one2many and not self.many2one:
			self.display_help("No files provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfastaFormatter.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description:Right now just converts b/n multi-line and one-line fasta formats, might add later")
		print("""
		-1,--one2many	: Path to fasta file to multi-line format
		-M,--many2one	: Path to fasta file to convert to one-line format
		-w,--width	: Characters per line for multi-line (default: 60)
		-o,--out	: Output file name (default=out.fas)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
