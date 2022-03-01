#!/usr/bin/python

import sys
import os
import getopt

import urllib.parse

def main():
	params = parseArgs()


	if params.gff:
		total_genes=0
		blacklist=list()
		current_gene=None
		current_sum=0
		for r in read_gff(params.gff):
			if r.type == "CDS":
				if params.gffid in r.attributes:
					name=r.attributes[params.gffid]
					if not current_gene:
						current_gene=name
					elif current_gene == name:
						current_sum += r.end-r.start
					else:
						#gene is done
						if current_sum%3 != 0:
							blacklist.append(current_gene)
						current_gene=None
						current_sum=0
						total_genes+=1
				else:
					print("ERROR: Given ID not present in GFF attributed field.")
					sys.exit()

		if len(blacklist)>0:
			print("Found", str(len(blacklist)), "of", str(total_genes),"genes not divisible by 3:")
			for g in blacklist:
				print(g)
		else:
			print("All",str(total_genes),"genes passed.")

#Function to split GFF attributes
def splitAttributes(a):
	ret = {}
	for thing in a.split(";"):
		stuff = thing.split("=")
		if len(stuff) != 2: continue #Eats error silently, YOLO
		key = stuff[0].upper()
		value = stuff[1].upper()
		ret[key] = value
	return ret

#Class for holding GFF Record data, no __slots__
class GFFRecord():
	def __init__(self, things):
		self.seqid = "NULL" if things[0] == "." else urllib.parse.unquote(things[0])
		self.source = "NULL" if things[1] == "." else urllib.parse.unquote(things[1])
		self.type = "NULL" if things[2] == "." else urllib.parse.unquote(things[2])
		self.start = "NULL" if things[3] == "." else int(things[3])
		self.end = "NULL" if things[4] == "." else int(things[4])
		self.score = "NULL" if things[5] == "." else float(things[5])
		self.strand = "NULL" if things[6] == "." else urllib.parse.unquote(things[6])
		self.phase = "NULL" if things[7] == "." else urllib.parse.unquote(things[7])
		self.attributes = {}
		if things[8] != "." and things[8] != "":
			self.attributes = splitAttributes(urllib.parse.unquote(things[8]))

	def getAlias(self):
		"""Returns value of alias if exists, and False if it doesn't exist"""
		if 'alias' in self.attributes:
			return self.attributes['alias']
		else:
			return False

#file format:
#1 per line:
#chr1:1-1000
#...
def read_regions(r):
	with open(r, 'w') as fh:
		try:
			for line in fh:
				line = line.strip() #strip leading/trailing whitespace
				if not line: #skip empty lines
					continue
				yield(ChromRegion(line))
		finally:
			fh.close()


#function to read a GFF file
#Generator function, yields individual elements
def read_gff(g):
	bad = 0 #tracker for if we have bad lines
	gf = open(g)
	try:
		with gf as file_object:
			for line in file_object:
				if line.startswith("#"): continue
				line = line.strip() #strip leading/trailing whitespace
				if not line: #skip empty lines
					continue
				things = line.split("\t") #split lines
				if len(things) != 9:
					if bad == 0:
						print("Warning: GFF file does not appear to be standard-compatible. See https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md")
						bad = 1
						continue
					elif bad == 1:
						sys.exit("Fatal error: GFF file does not appear to be standard-compatible. See https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md")
				#line = utils.removeURL(line) #Sanitize any URLs out
				rec = GFFRecord(things)

				yield(rec)
	finally:
		gf.close()

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
                    options, remainder = getopt.getopt(sys.argv[1:], 'hg:o:i:p:r:', \
			["help", "gff=", "out=", "printAll", "remove", "id="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.gff=None
		self.out=None
		self.gffid="protein_id"
		self.remove=False
		self.printAll=False


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
			elif opt=="gff" or opt=="g":
				self.gff=arg
			elif opt=="o" or opt=="out":
				self.out=arg
			elif opt=="r" or opt=="remove":
				self.remove=True
			elif opt=="p" or opt=="printAll":
				self.printAll=True
			elif opt=="i" or opt=="id":
				self.gffid=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.gff:
			self.display_help("Must provide a GFF file")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ncheck_cds.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description:Checks that the sum of CDS blocks for a gene is divisible by 3")
		print("""
		-g,--gff	: GFF File (may contain non-CDS regions but these will be skipped)
		-o,--out	: Output file for new CDS (default=None)
		-r,--remove	: Remove genes not divisible by 3
		-i,--id	: Field in GFF attributes (last column) giving the identifier to use
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
