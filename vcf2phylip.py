#!/usr/bin/python

import re
import sys
import os
import vcf
import getopt

def main():
	params = parseArgs()

	data = dict()

	if params.vcf:
		#for each record in VCF
		for record in read_vcf(params.vcf):
			for call in record.samples:
				#Get consensus base call
				cons = None
				if call.gt_bases:
					l = (call.gt_bases).split("/")
					cons = reverse_iupac(listToSortUniqueString(l))
				else:
					cons = "N"
				if cons:
					if call.sample in data:
						data[call.sample].append(cons)
					else:
						data[call.sample] = list()
						data[call.sample].append(cons)
				else:
					print ("Uh oh! No consensus called for %s, something is wrong"%call)

		#Print dict to phylip file
		with open(params.out, 'w') as fh:
			try:

			except IOError:
				print("Could not write to file ",params.out)
				sys.exit(1)
			finally:
				fh.close()

	else:
		print("Error: No VCF file provided")
		sys.exit(1)


#Read VCF variant calls
#Generator function, yields each locus
def read_vcf(v):

	try:
		vfh = vcf.Reader(filename=v)
	except IOError as err:
		print("I/O error({0}): {1}".format(err.errno, err.strerror))
	except:
		print("Unexpected error:", sys.exec_info()[0])

	chrom = ""
	recs = []
	added = 0
	for rec in vfh:
		if not rec.FILTER:
			yield(rec)

#Function to return sorted unique string from list of chars
def listToSortUniqueString(l):
	sl = sorted(set(l))
	return(str(''.join(sl)))

#Function to translate a string of bases to an iupac ambiguity code
def reverse_iupac(char):
	char = char.upper()
	if "-" in char:
		return("-")
	elif "N" in char:
		return("N")
	elif "." in char:
		return(".")
	else:
		iupac = {
			'A':'A',
			'N':'N',
			'-':'-',
			'C':'C',
			'G':'G',
			'T':'T',
			'AG':'R',
			'CT':'Y',
			'AC':'M',
			'GT':'K',
			'AT':'W',
			'CG':'S',
			'CGT':'B',
			'AGT':'D',
			'ACT':'H',
			'ACG':'V',
			'ACGT':'N'
		}
		return iupac[char]

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'v:o:h', \
			["vcf=","help","out="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.vcf=None
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
			if opt in ('v', 'vcf'):
				self.vcf = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.vcf:
			self.display_help("\nError: Missing required input file <-v,--vcf>")


		if self.out:
			self.out = self.out + ".phy"
		else:
			self.out = "out.phy"


	def display_help(self, message=None):
		if message is not None:
			print (message)
		print ("\nvcf2phylip.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-v </path/to/vcf>\n")
		print ("Description: Extract SNPs from a VCF file and outputs as concatenated Phylip")

		print("""
	Arguments:
		-v,--vcf	: VCF input file
		-o,--out	: Prefix for output file <default = ./out>
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
