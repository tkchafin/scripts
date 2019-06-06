#!/usr/bin/python

import sys
import os
import getopt
import io
import glob
from Bio import AlignIO
from Bio import SeqIO


def main():
	params = parseArgs()

	if params.files:
		alignments = list()
		lookup = dict()
		parts = dict()
		for f in glob.iglob(params.files):
			#Read file into dict of dicts
			this_one = dict()
			this_length = None

			#for each tuple in NEXUS file
			print("Reading file:",f)

			#Read nexus as DNAalignment
			this_one = DNAalignment(f, "nexus")
			if this_one.length < params.min:
				print("Alignment",this_one.name,"has length",this_one.length,", skipping.")
				continue


			# if params.part:
			# 	this_one.loadPartitions(f)

			#add names to lookup table
			for name in this_one.alignment.keys():
				lookup[name] = "NULL"

			#Add this alignment to our alignment list
			alignments.append(this_one)

		#For each alignment, print individual NEXUS if required
		if params.ind_files:
			for aln in alignments:
				#print NEXUS header information
				new_nex = aln.name + "_filled.nex"
				print("Writing",new_nex)
				print(new_nex)
				with open(new_nex, 'w') as fh:
					try:
						slen = aln.length
						header = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=" + str(len(lookup.keys())) + " NCHAR=" + str(slen) + ";\n"
						header = header + "FORMAT DATATYPE=DNA MISSING=? GAP=-;\n\nMATRIX\n"
						fh.write(header)
						for samp in sorted(lookup.keys()):
							#if sample not in alignment, write Ns
							if samp not in aln.alignment.keys():
								l = str(samp) + "\t" + str(Nrepeats("?", aln.length)) + "\n"
								fh.write(l)
							else:
								l = str(samp) + "\t" + str(aln.alignment[samp]) + "\n"
								fh.write(l)
						last = ";\nEND;\n"
						fh.write(last)
					except IOError as e:
						print("Could not write to file:",new_nex)
						print("Caught IOError: ",e)
						sys.exit(1)
					except Exception as e:
						print("Unexpected error writing to file:",new_nex)
						print("Caught Exception: ",e)
						sys.exit(1)
					finally:
						fh.close()

		#write concatenated file with CHARSETS
		concat = "concat.nex"
		with open(concat, 'w') as fh:
			try:
				#write header
				print ("Writing concatenated alignment to concat.nex...")
				full_len = 0
				for a in alignments:
					full_len += a.length
				header = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=" + str(len(lookup.keys())) + " NCHAR=" + str(full_len) + ";\n"
				header = header + "FORMAT DATATYPE=DNA MISSING=? GAP=-;\n\nMATRIX\n"
				fh.write(header)

				for samp in sorted(lookup.keys()):
					l = samp + "\t"
					for aln in alignments:
						#if sample not in alignment, write Ns
						if samp not in aln.alignment.keys():
							l += str(Nrepeats("?", aln.length))
						else:
							l += str(aln.alignment[samp])
					l += "\n"
					fh.write(l)
				l = "\n;END;\n"
				fh.write(l)

				#write charsets
				print("Adding CHARSET blocks...")
				out = "\nBEGIN sets;\n"
				cur = 0
				for a in alignments:
					out += "\t" + "CHARSET " + str(a.name) + " = " + str(cur + 1) + "-" + str(cur + a.length) + ";\n"
					cur += a.length
				out += "END;\n\n"
				fh.write(out)

			except IOError as e:
				print("Could not write to file:",concat)
				print("Caught IOError: ",e)
				sys.exit(1)
			except Exception as e:
				print("Unexpected error writing to file:",concat)
				print("Caught Exception: ",e)
				sys.exit(1)
			finally:
				fh.close()


		print("\nDone!\n")

	else:
		print("No inputs provided.\n")
		sys.exit(0)


def Nrepeats(pattern, N):
	ret = ""
	for i in range(int(N)):
		ret = ret + str(pattern)
	return(ret)



#Object to hold an alignment
class DNAalignment():
	def __init__(self, aln_file, aln_type):
		self.length = int()
		self.alignment = dict()
		self.filepath = aln_file
		self.name = self.getName()
		if aln_type.lower() in ["nexus", "nex"]:
			for sample in self.readNexus():
				#Add to this alignment
				if sample[0] not in self.alignment.keys():
					self.alignment[sample[0]] = sample[1]
					if not self.length:
						self.length = len(sample[1])
					else:
						if len(sample[1]) != self.length:
							print("Sequences not of equal length:",f)
							sys.exit(1)
	def getName(self):
		return(os.path.splitext(os.path.basename(self.filepath))[0])

	#method to read alignment from NEXUS
	def readNexus(self):
		with open(self.filepath, 'r') as fh:
			try:
				start = False
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#line = line.replace(" ","")

					if line.lower() == "matrix":
						start = True
						continue

					if start: #we're in the matrix!
						if line in [";", "END;", "end;"]:
							start = False
							break
						else:
							line = line.replace("\'","")
							line = line.replace("\"","")
							stuff = line.split()
							yield([stuff[0],stuff[1]])
			except IOError as e:
				print("Could not read file:",e)
				sys.exit(1)
			except Exception as e:
				print("Unexpected error:",e)
				sys.exit(1)
			finally:
				fh.close()

	#Function to read partitions from a NEXUS-formatted SETS block
	# def loadPartitions(self, f):
	# 	with open(f, 'r') as fh:
	# 		try:
	# 			start = False
	# 			for line in fh:
	# 				line = line.strip()
	# 				if not line:
	# 					continue
	# 				#line = line.replace(" ","")
	#
	# 				if "begin" and "sets" in line.lower():
	# 					start = True
	# 					continue
	#
	# 				if start: #we're in the charsets block
	# 					if "charset" in line.lower():
	#
	#
	# 		except IOError as e:
	# 			print("Could not read file:",e)
	# 			sys.exit(1)
	# 		except Exception as e:
	# 			print("Unexpected error:",e)
	# 			sys.exit(1)
	# 		finally:
	# 			fh.close()

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hf:ipmM:', \
			["help", "files=","ind", "part", "miss",'min'])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.files = None
		self.ind_files = False
		self.part = False
		self.recode = False
		self.min = 1


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
			if opt in ('f','files'):
				self.files = arg
			elif opt in ('i','ind'):
				self.ind_files = True
			elif opt in ('p','part'):
				self.part = True
			elif opt in ('m', 'miss'):
				self.recode = True
			elif opt in ('M', 'min'):
				self.min=int(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.files:
			self.display_help("No files provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nconcatenate.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Concatenates NEXUS gene alignments and fills in gaps for missing taxa")
		print("""
	Options:
		-f,--files	: \"/path/to/*.nex\" MUST BE QUOTED, NO \'~\'!!
		-i,--ind	: Toggle on to print individual NEXUS files [default=off]
		-p,--part	: Toggle on to retain existing partition information (CHARSETS)
		-m,--miss	: Toggle on to code missing data as "?"
			NOTE: This will also convert terminal gaps as "?"
		-M,--min	: Minimum alignment length to keep an alignment
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
