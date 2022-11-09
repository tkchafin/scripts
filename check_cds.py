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
		gene_sums=dict()

		if params.correct:
			print("WARNING: You have chosen to force non-conforming CDS to be divisible by 3. Note that this assumes CDS records are sorted")
		if params.remove:
			print("WARNING: Removing all non-conforming CDS from output file")

		for r in read_gff(params.gff):
			if r.type == "CDS":
				if params.gffid in r.attributes:
					name=r.attributes[params.gffid]
					if not current_gene:
						current_gene=name
						current_sum += (r.end-r.start)+1
					elif current_gene == name:
						current_sum += (r.end-r.start)+1
					else:
						#gene is done
						if current_sum%3 != 0:
							blacklist.append(current_gene)
						gene_sums[current_gene]=current_sum
						current_gene=name
						current_sum=0
						total_genes+=1
				else:
					print("ERROR: Given ID not present in GFF attributed field.")
					sys.exit()
		if current_sum>0:
			if current_sum%3 != 0:
				blacklist.append(current_gene)
			gene_sums[current_gene]=current_sum
			current_gene=name
			current_sum=0
			total_genes+=1

		if len(blacklist)>0:
			print("Found", str(len(blacklist)), "of", str(total_genes),"genes not divisible by 3:")
			for g in blacklist:
				print(g)

			if params.out:
				print("Writing new output GFF (CDS records only):",params.out)
				with open(params.out, 'w') as fh:
					last_line=None
					current_gene=None
					if len(blacklist)>0:
						lookup = dict(zip(blacklist, blacklist))
						for r in read_gff(params.gff):
							if params.gffid in r.attributes:
								name=r.attributes[params.gffid]

								if r.type != "CDS":
									if params.printAll:
										if params.remove:
											if name not in blacklist:
												fh.write(r.as_string())
												fh.write("\n")
										else:
											fh.write(r.as_string())
											fh.write("\n")
								else:
									#if continuing gene, check if blacklisted
									if not current_gene:
										current_gene=name

									if current_gene == name:
										#print(name)
										if name in blacklist:
											#print("bad")
											#skip if params.remove
											if params.correct:
												if last_line:
													fh.write(last_line)
													fh.write("\n")
													last_line=r.as_string()
											last_line=r.as_string()
										else:
											fh.write(r.as_string())
											fh.write("\n")
									else:
										if current_gene in blacklist:
											#print(name)
											#print("bad")
											if params.correct:
												#print("correcting",current_gene)
												if last_line:
													#print("last_line:")
													#print(last_line)
													corrected=GFFRecord(last_line.split("\t"))
													#print(corrected.as_string())
													#print("GENE: ", corrected.attributes[params.gffid])
													#print("GENE SUM: ",gene_sums[corrected.attributes[params.gffid]])
													#print("DIFFERENCE:", gene_sums[corrected.attributes[params.gffid]]%3)
													diff=gene_sums[corrected.attributes[params.gffid]]%3
													corrected.end-=diff
													#print(corrected.as_string())
													fh.write(r.as_string())
													fh.write("\n")
										current_gene=name

						if current_gene in blacklist:
							#print(name)
							#print("bad")
							if params.correct:
								#print("correcting",current_gene)
								if last_line:
									#print("last_line:")
									#print(last_line)
									corrected=GFFRecord(last_line.split("\t"))
									#print(corrected.as_string())
									#print("GENE: ", corrected.attributes[params.gffid])
									#print("GENE SUM: ",gene_sums[corrected.attributes[params.gffid]])
									#print("DIFFERENCE:", gene_sums[corrected.attributes[params.gffid]]%3)
									diff=gene_sums[corrected.attributes[params.gffid]]%3
									corrected.end-=diff
									#print(corrected.as_string())
									fh.write(r.as_string())
									fh.write("\n")

		else:
			print("All",str(total_genes),"genes passed.")

#Function to split GFF attributes
def splitAttributes(a):
	ret = {}
	for thing in a.split(";"):
		stuff = thing.split("=")
		if len(stuff) != 2: continue #Eats error silently, YOLO
		key = stuff[0]
		value = stuff[1]
		ret[key] = value
	return ret

#Class for holding GFF Record data, no __slots__
class GFFRecord():
	def __init__(self, things):
		self.seqid = "." if things[0] == "." else urllib.parse.unquote(things[0])
		self.source = "." if things[1] == "." else urllib.parse.unquote(things[1])
		self.type = "." if things[2] == "." else urllib.parse.unquote(things[2])
		self.start = "." if things[3] == "." else int(things[3])
		self.end = "." if things[4] == "." else int(things[4])
		self.score = "." if things[5] == "." else float(things[5])
		self.strand = "." if things[6] == "." else urllib.parse.unquote(things[6])
		self.phase = "." if things[7] == "." else urllib.parse.unquote(things[7])
		self.attributes = {}
		if things[8] != "." and things[8] != "":
			self.attributes = splitAttributes(urllib.parse.unquote(things[8]))

	def getAlias(self):
		"""Returns value of alias if exists, and False if it doesn't exist"""
		if 'alias' in self.attributes:
			return self.attributes['alias']
		else:
			return False

	def as_string(self):
		record=str(self.seqid) + "\t" + str(self.source) + "\t" + str(self.type)
		record=record + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + str(self.score)
		record=record + "\t" + str(self.strand) + "\t" + str(self.phase) + "\t"
		first=True
		for key in self.attributes.keys():
			if first:
				record=record + str(key) + "=" + str(self.attributes[key])
				first=False
			else:
				record=record + ";" + str(key) + "=" + str(self.attributes[key])
		return(record)


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
                    options, remainder = getopt.getopt(sys.argv[1:], 'hg:o:i:prc', \
			["help", "gff=", "out=", "printAll", "remove", "id=", "correct"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.gff=None
		self.out=None
		self.gffid="protein_id"
		self.remove=False
		self.correct=False
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
			elif opt=="p" or opt=="printAll":
				self.printAll=True
			elif opt=="i" or opt=="id":
				self.gffid=arg
			elif opt=="c" or opt=="correct":
				self.correct=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.gff:
			self.display_help("Must provide a GFF file")
		if not self.correct:
			self.remove=True
		if self.correct and not self.out:
			self.display_help("Cannot use --correct without -o/--out")

	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ncheck_cds.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Checks that the sum of CDS blocks for a gene is divisible by 3")
		print("""
		-g,--gff	: GFF File (may contain non-CDS regions but these will be skipped)
		-o,--out	: Output file for new CDS (default=None)
		-c,--correct	: Walk back end coordinate for genes not summing to %3 (only relevant if --out)
		-i,--id	: Field in GFF attributes (last column) giving the identifier to use
		-p,--printAll:	Print all fields (not just CDS)(only relevant if --out)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
