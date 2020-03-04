#!/usr/bin/python

import sys
import os
import getopt

def main():
	params = parseArgs()
	locnum=1
	skipped=0
	contents=""
	print("Minimum allowable alignment length:",params.minlen)
	for file in os.listdir(params.fasdir):
		if file.endswith(".fas") or file.endswith(".fasta") or file.endswith(".fsa"):
			aln=dict()
			tax=0
			aln_len=0
			skip=False
			for line in read_fasta(params.fasdir + "/" + file):
				aln[line[0]] = line[1].replace("-", "N")
				tax=tax+1
				aln_len=len(line[1])
				if aln_len < params.minlen:
					skip=True
					continue
			if skip:
				skipped=skipped+1
				continue
			contents=contents+"locus"+str(locnum)+" "+str(tax)+" "+str(aln_len)+"\n"
			locnum=locnum+1
			#print(locnum)
			for samp in sorted(aln):
				#print(samp)
				contents=contents+str(samp)+" "+aln[samp]+"\n"
			contents=contents+"\n"
	#print(contents)

	print("Skipped alignments smaller than minimum length:",skipped)
	print("Total alignments passing filtering:",locnum)
	
	ofh=open(params.out, "w")
	header=str(locnum)+"\n"
	ofh.write(header)
	ofh.write(contents)
	ofh.close()

#Read genome as FASTA. FASTA header will be used
#This is a generator function
#Doesn't matter if sequences are interleaved or not.
def read_fasta(fas):
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


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hf:o:m:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.fasdir=None
		self.out="gphocs_input.txt"
		self.minlen=500

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
			elif opt == "f":
				self.fasdir=arg
			elif opt=="o":
				self.out=arg
			elif opt=="m":
				self.minlen=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.fasdir:
			self.display_help("No files provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfasta2gphocs.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Converts a set of separate FASTA-formatted gene alignments to g-phocs sequence file format")
		print("""
	Arguments:
	-f	: Directory containing FASTA files
	-o	: Output file name
	-m	: Minimum alignment length (default=500)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
