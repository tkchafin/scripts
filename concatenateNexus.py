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
		
		#if -t, read and validate popmap
		if params.popmap is not None:
			popmap=parsePopmap(params.popmap)
			for ind in lookup.keys():
				if ind not in popmap:
					print("Sample", ind, "not found in popmap. Treating as separate pop.")
					popmap[ind] = ind
			blacklist=list()
			for ind in popmap.keys():
				if ind not in lookup:
					print("Sample", ind, "found in popmap but not in data. Deleting it.")
					blacklist.append(ind)
			for ind in blacklist:
				del popmap[ind]
			#make flattened popmap
			flatmap=make2Dpopmap(popmap)
		
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
						if params.popmap is not None:
							for pop in sorted(list(flatmap.keys())):
								for samp in flatmap[pop]:
									l = samp + "\t"
									for aln in alignments:
										#if sample not in alignment, write Ns
										if samp not in aln.alignment.keys():
											l += str(Nrepeats("?", aln.length))
										else:
											l += str(aln.alignment[samp])
									l += "\n"
									fh.write(l)
						else:
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
						

						partline1 = "\nbegin sets;"
						partline2 = "\ttaxpartition popmap ="
						fh.write("{}\n{}\n".format(partline1, partline2))
						
						#if -t, write taxpartition block 
						if params.popmap is not None:
							popnum = 1
							sample_number = 1
							range_begin = 1
							prev = None
							for pop in sorted(list(flatmap.keys())):
								start = sample_number
								for ind in flatmap[pop]:
									sample_number += 1
								end=sample_number
								if prev is not None:
									write_taxpart(prev[0], fh, tuple([prev[1], prev[2]-1]))
								prev=[pop, start, end]
							fh.write("\t\t{}\t:	{}-{};\n".format(str(prev[0]), str(prev[1]), str(prev[2]-1)))
							fh.write("end;\n")
							
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
				
				if params.popmap is not None:
					for pop in sorted(list(flatmap.keys())):
						for samp in flatmap[pop]:
							l = samp + "\t"
							for aln in alignments:
								#if sample not in alignment, write Ns
								if samp not in aln.alignment.keys():
									l += str(Nrepeats("?", aln.length))
								else:
									l += str(aln.alignment[samp])
							l += "\n"
							fh.write(l)
							
				else:
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
				
				#end of data block
				l = "\n;END;\n"
				fh.write(l)
				
				partline1 = "\nbegin sets;"
				partline2 = "\ttaxpartition popmap ="
				fh.write("{}\n{}\n".format(partline1, partline2))

				#write charsets
				print("Adding CHARSET blocks...")
				out = "\nBEGIN sets;\n"
				cur = 0
				for a in alignments:
					out += "\t" + "CHARSET " + str(a.name) + " = " + str(cur + 1) + "-" + str(cur + a.length) + ";\n"
					cur += a.length
				out += "END;\n\n"
				fh.write(out)
				
				#if -t, write taxpartition block 
				if params.popmap is not None:
					popnum = 1
					sample_number = 1
					range_begin = 1
					prev = None
					for pop in sorted(list(flatmap.keys())):
						start = sample_number
						for ind in flatmap[pop]:
							sample_number += 1
						end=sample_number
						if prev is not None:
							write_taxpart(prev[0], fh, tuple([prev[1], prev[2]-1]))
						prev=[pop, start, end]
					fh.write("\t\t{}\t:	{}-{};\n".format(str(prev[0]), str(prev[1]), str(prev[2]-1)))
					fh.write("end;\n")

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

def write_taxpart(pattern, outfile, range_tupl):
	outfile.write("\t\t{}\t:    {}-{},\n".format(str(pattern), str(range_tupl[0]), str(range_tupl[1])))


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

#Makes a dict of lists from a popmap
def make2Dpopmap(p):
	ret = dict()
	for s in p:
		if p[s] not in ret:
			ret[p[s]] = list()
		ret[p[s]].append(s)
	return(ret)

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
			options, remainder = getopt.getopt(sys.argv[1:], 'hf:ipmM:t:', \
			["help", "files=","ind", "part", "miss",'min', "tax="])
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
		self.popmap=None


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
			elif opt == "t" or opt == "tax":
				self.popmap=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.files:
			self.display_help("No files provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nconcatenateNexus.py\n")
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
		-t, --tax	: Popmap to write taxpartition
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
	main()
