#!/usr/bin/python

import sys
import os
import getopt

def main():
	params = parseArgs()
	
	f = open(params.out, 'w')
	
	with open(params.vcf, "r") as vcf:
		for line in vcf:
			line=line.strip()
			#directly transfer header lines
			if line[0] == "#":
				if line[1] != "#":
					f.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n")
				f.write(line)
				f.write("\n")
				
			else:
				fields=line.split("\t")
				ref=get_index(fields[3].split(","))
				alt=get_index(fields[4].split(","))
				#if biallelic filter on and site has >2 alleles, skip
				if params.biallelic == True:
					if (len(ref) + len(alt) > 2):
						continue
				if fields[8]=="GT:DP:CATG":
					fields[8]="GT:DP:AD"
				else:
					print("Something wrong with VCF. Field 8 should be GT:DP:CATG")
					sys.exit()
				for idx, sample in enumerate(fields[9:]):
					fixed=str(fix_sample(sample, ref, alt))+":"
					fields[idx+9] = fixed
					#print(sample, " -- ", fixed)
				f.write("\t".join(fields))
				f.write("\n")
		vcf.close()
	f.close()

def fix_sample(sample, ref, alt):
	fields=sample.split(":")
	catg=fields[2].split(",")
	ad=list()
	for r in ref:
		ad.append(catg[r])
	for a in alt:
		ad.append(catg[a])
	fields[2]=",".join(ad)
	return(":".join(fields))
	
def get_index(char):
	ret=list()
	for c in char:
		if c.lower()=="c":
			ret.append(0)
		elif c.lower()=="a":
			ret.append(1)
		elif c.lower()=="t":
			ret.append(2)
		elif c.lower()=="g":
			ret.append(3)
		else:
			print("Unrecognized character",char)
			sys.exit()
	return(ret)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hv:o:b', \
			["help", "vcf=", "out=", "biallelic"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
		self.out="polyrad.vcf"
		self.biallelic=False


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
			elif opt=="vcf" or opt=="v":
				self.vcf=arg
			elif opt=="out" or opt=="o":
				self.out=arg
			elif opt=="biallelic" or opt=="b":
				self.biallelic=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("Need an ipyrad VCF file")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nipyrad2polyrad.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description:Converts the ipyrad VCF to a format usable for polyRAD")
		print("""
		-v,--vcf	: VCF input with ipyrad "CATG" field
		-b,--biallelic	: [Boolean] Toggle to skip non-biallelic sites
		-o,--out	: Output file name (default=polyrad.vcf)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
