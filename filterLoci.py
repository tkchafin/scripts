#!/usr/bin/python

import re
import sys
import os
import getopt
import random


def main():
	params = parseArgs()
	
	whitelist=list()
	if params.whitelist:
		with open(params.whitelist, "r") as wl:
			for line in wl:
				whitelist.append(int(line.strip()))
	if params.input:
		if params.popmap:
			print("Parsing popmap file...")
			pop_assign = parsePopmap(params.popmap)

		allLoci = list()
		print("Requiring a minimum of", params.samples, "per locus")
		locNum = 0
		with open(params.input, "r") as fh:
			try:
				aln_d = dict()
				numSamp = int(0)
				numPIS = int(0)

				passed=0
				sampFail=0
				varsFail=0
				snaqFail=0
				listFail=0


				if params.snaqr:
					print("USING SNAQR PRESET. ONLY LOCI WITH AT LEAST 1 SAMPLE PER POP WILL BE KEPT")
				if params.snaq or params.snaqm:
					bestInds = getBestIndLoci(params.input, pop_assign, params.pis)
					print("USING SNAQ PRESET. THE FOLLOWING INDIVIDUALS WILL BE SAMPLED AS REPRESENTATIVES:")
					for pop in bestInds:
						print("%s: %s"%(pop, bestInds[pop]))
				if params.pomo:
					print("USING POMO PRESET. ONLY LOCI WITH AT LEAST 1 SAMPLE PER POP WILL BE KEPT")

				passing=dict()

				for line in fh:
					line = line.strip()
					if not line:
						continue
					if line and \
					(line[0].isalpha() or \
					line[0].isdigit() or \
					line[0] != "/"):
						allLoci.append(line)
						stuff = line.split()
						if line[0] != "/":
							identifier = stuff[0].replace(">", "")
						else:
							identifier = stuff[0]
						seqs = stuff[1]
						aln_d[identifier] = seqs
						#print(len(seqs))
						numSamp += 1
					elif line[0] == "/": #alignment end. grab numPIS and parse locus
						#print("loc", locNum)
						numPIS = line.count("*")
						allLoci.append(line)
						locNum += 1
						
						if len(whitelist) > 0:
							if locNum not in whitelist:
								listFail+=1
								aln_d=dict()
								numSamp=0
								numPIS=0
								continue
						#print(numSamp)
						if numPIS >= params.pis and numSamp >= params.samples:

							passed += 1
							#Write to desired location
							if params.fas:
								#f = str(params.out) + "_" + str(locNum) + ".fasta"
								#writeFasta(aln_d, f)
								passing[locNum]=aln_d
								aln_d=dict()
								numPIS=0
								numSamp=0

							if params.nex:
								n = str(params.out) + "_" + str(locNum) + ".nex"
								if params.snaq:
									aln_new = sampledSnaq(aln_d,bestInds)
									if aln_new:
										#dict2nexus(n, aln_new)
										passing[locNum]=aln_new
										aln_d=dict()
										numPIS=0
										numSamp=0
									else:
										snaqFail+=1
										passed -= 1
										aln_d=dict()
										numPIS=0
										numSamp=0
										continue
								elif params.snaqm:
									aln_new = sampledSnaqM(aln_d,bestInds)
									if len(aln_new) >= params.samples:
										#dict2nexus(n, aln_new)
										passing[locNum]=aln_new
										aln_d=dict()
										numPIS=0
										numSamp=0
									else:
										snaqFail+=1
										passed -= 1
										numSamp=0
										numPIS=0
										aln_d=dict()
										continue
								elif params.snaqr:
									aln_new = randomSnaq(aln_d,pop_assign)
									if aln_new:
										#dict2nexus(n, aln_new)
										passing[locNum]=aln_new
										aln_d=dict()
										numPIS=0
										numSamp=0
									else:
										snaqFail+=1
										passed -= 1
										numSamp=0
										numPIS=0
										aln_d=dict()
										continue
								elif params.pomo:
									aln_new = sampledPop(aln_d, pop_assign, params.pomo)
									if aln_new:
										#dict2nexus(n, aln_new)
										passing[locNum]=aln_new
										aln_d=dict()
										numPIS=0
										numSamp=0
									else:
										snaqFail+=1
										passed -= 1
										numSamp=0
										numPIS=0
										aln_d=dict()
										continue
								else:
									passing[locNum]=aln_d
									dict2nexus(n, aln_d)
									aln_d=dict()
									numPIS=0
									numSamp=0
							numSamp=0
							numPIS=0
						else:
							if numPIS < params.pis:
								varsFail += 1
								aln_d=dict()
								numPIS=0
								numSamp=0
								continue
							if numSamp < params.samples:
								sampFail += 1
								aln_d=dict()
								numPIS=0
								numSamp=0
								continue
						aln_d=dict() #clear alignment
						numSamp=0
						numPIS=0

			except IOError as e:
				print("Couldn't read file %s: %s" %(params.input,e))
				sys.exit(1)
			except Exception as e:
				print("Unexpected error reading file %s: %s" %(params.input, e))
				sys.exit(1)
			finally:
				fh.close()

		# if params.loci:
		#   l = str(params.out) + "_filtered.loci"
		#   writeStuff(allLoci, l)
		print("Number of loci passing filters: %s (of %s total)" %(passed, locNum))
		print("\t",sampFail,"loci failed for too few individuals.")
		print("\t",varsFail,"loci failed for too few parsimony-informative sites\n")
		if params.snaq or params.snaqr:
			print("\t",snaqFail,"loci failed SNAQ-sampling\n")
		
		if params.rand:
			print("Randomly sampling",params.rand,"passing loci...")
			stuff = list(passing.keys())
			random.shuffle(stuff)
			max = params.rand
			new = dict()
			count = 0
			for key in stuff:
				if count >= max:
					break
				new[key] = passing[key]
				count += 1
			passing = new

		if params.fas:
			for aln in passing:
				f = str(params.out) + "_" + str(aln) + ".fasta"
				writeFasta(passing[aln], f)
		if params.nex:
			for aln in passing:
				pass
				#print("Writing",n)
				n = str(params.out) + "_" + str(aln) + ".nex"
				dict2nexus(n, passing[aln])
				

	else:
		print("No input provided.")
		sys.exit(1)

def sampledSnaqM(aln, best):
	ret = dict()
	for pop in best.values():
		if pop in aln:
			ret[pop] = aln[pop]
	return(ret)

def sampledSnaq(aln, best):
	ret = dict()
	for pop in best.values():
		if pop in aln:
			ret[pop] = aln[pop]
		else:
			return(False)
	return(ret)

def randomSnaq(aln,popmap):
	pop_enum = set(popmap.values())
	#print(pop_enum)
	ret=dict()

	for pop in pop_enum:
		options = list()
		for samp in aln.keys():
			if samp in popmap.keys() and pop == popmap[samp]:
				options.append(samp)
		if len(options) < 1:
			return(False)
		elif len(options) == 1:
			ret[pop] = aln[options[0]]
		else:
			ret[pop] = aln[random.choice(options)]
	return(ret)

def sampledPop(aln,popmap, num):
	pop_enum = set(popmap.values())
	#print(pop_enum)
	ret=dict()

	for pop in pop_enum:
		options = list()
		for samp in aln.keys():
			if samp in popmap.keys() and pop == popmap[samp]:
				options.append(samp)
		if len(options) < 1:
			return(False)
		elif len(options) < num:
			for opt in options:
				ret[pop] = aln[opt]
		else:
			random.shuffle(options)
			for opt in options[0:num-1]:
				ret[pop] = aln[opt]
	return(ret)


#Function returns dict of highest coverage individual for each pop in pop_assign
def getBestIndLoci(infile, popmap, thresh):
	pop_enum = set(popmap.values())

	indcovs=dict()
	for ind in popmap.keys():
		indcovs[ind] = 0

	#Get individual coverage by scanning loci file
	with open(infile, "r") as fh:
		try:
			seen = list()
			num=0
			for line in fh:
				line = line.strip()
				if not line:
					continue
				if line[0] == ">" or \
				line[0].isalpha() or \
				line[0].isdigit():
					if line[0] == ">":
						line = line.replace(">","")
					stuff = line.split()
					identifier = stuff[0].replace(" ", "")
					seen.append(identifier)

				elif line[0] == "/": #alignment end. grab numPIS and parse locus
					num+=1
					numPIS = line.count("*")
					if numPIS >= thresh: #locus passes
						for samp in seen:
							if samp in indcovs:
								indcovs[samp] += 1
							else:
								indcovs[samp] = 1
					seen = list()
		except IOError as e:
			print("Couldn't read file %s: %s" %(infile,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s" %(infile, e))
			sys.exit(1)
		finally:
			fh.close()

	#Figure out which samples per pop have the highest coverage
	bestDict = dict()
	for pop in pop_enum:
		best = None
		bestcount = None
		for ind in popmap.keys():
			if popmap[ind] == pop:
				if best:
					if indcovs[ind] > bestcount:
						best = ind
						bestcount = indcovs[ind]
				else:
					best = ind
					bestcount = indcovs[ind]
		bestDict[pop] = best

	return(bestDict)

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

#function writes everything in list to a file
def writeStuff(stuff, f):
	with open(f, "w") as fh:
		try:
			for line in stuff:
				out = str(line) + "\n"
				fh.write(out)
		except IOError as e:
			print("Couldn't write to file %s: %s" %(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error writing to file %s: %s" %(f, e))
			sys.exit(1)
		finally:
			fh.close()

#function writes sequences in dictionary d to fasta-file f
def writeFasta(d, f):
	with open(f, "w") as ffh:
		try:
			for samp, seq in d.items():
				l = ">" + str(samp) + "\n" + seq + "\n"
				ffh.write(l)
		except IOError as e:
			print("Couldn't write to file %s: %s" %(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error writing to file %s: %s" %(f, e))
			sys.exit(1)
		finally:
			ffh.close()

#Function to write an alignment as DICT to NEXUS
def dict2nexus(nex, aln):
	with open(nex, 'w') as fh:
		try:
			slen = getSeqLen(aln)
			header = "#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX=" + str(len(aln)) + " NCHAR=" + str(slen) + ";\n"
			header = header + "FORMAT DATATYPE=DNA MISSING=N GAP=-;\n\nMATRIX\n"
			fh.write(header)
			for seq in aln:
				sline = str(seq) + "\t" + aln[seq] + "\n"
				fh.write(sline)
			last = ";\nEND;\n"
			fh.write(last)
		except IOError:
			print("Could not read file ",nex)
			sys.exit(1)
		finally:
			fh.close()

#Goes through a dict of sequences and get the alignment length
#returns minimum length
def getSeqLen(aln):
	length = None
	for key, val in aln.items():
		if not length:
			length = len(val)
		else:
			if length != len(val):
				print("Warning! Sequences in alignment not of equal length! Writing to BAD.fasta")
				writeFasta(aln, "BAD.fasta")
				#sys.exit(0)
	return(length)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'i:s:p:ho:nflP:X:SRMV:r:w:', \
			["input=", "samples=", "pis=","help","out=", "nex", "fas", "loci",
			"popmap=", "snaq", "exclude=","snaqR","snaqr", 'snaqM', 'snaqm', 'pomo', 'concatFas', "rand="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.input=None
		self.popmap=None
		self.samples=1
		self.pis=1
		self.out="loc"
		self.whitelist=None

		self.nex=False
		self.fas=False
		self.concatFas=False
		self.loci=False
		self.rand=None

		self.snaq=False
		self.snaqr=False
		self.snaqm=False
		self.pomo=None
		self.exclude=list()


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
			if opt in ('i', 'input'):
				self.input = arg
			elif opt in ('p', 'pis'):
				self.pis=int(arg)
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ('s', 'samples'):
				self.samples = int(arg)
			elif opt in ('n', 'nex'):
				self.nex=True
			elif opt in ('f', 'fas'):
				self.fas=True
			elif opt == 'l' or opt == 'loci':
				self.loci=True
			elif opt == 'S' or opt == 'snaq':
				self.snaq=True
			elif opt =='r' or opt == 'rand':
				self.rand=int(arg)
			elif opt == 'R' or opt == 'snaqr' or opt == 'snaqR':
				self.snaqr=True
			elif opt == 'M' or opt == 'snaqm' or opt == 'snaqM':
				self.snaqm=True
			elif opt == 'X' or opt == 'exclude':
				self.exclude=arg.split(",")
			elif opt in ('P','popmap'):
				self.popmap=arg
			elif opt in ('V', 'pomo'):
				self.pomo=int(arg)
			elif opt == "concatFas":
				self.concatFas=True
			elif opt=="w":
				self.whitelist=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.input:
			self.display_help(".loci file must be provided (-l, --loci)")

		if self.snaq or self.snaqr or self.pomo:
			self.nex = True
			if not self.popmap:
				self.display_help("You must provide a popmap with presets")
		if self.snaq and self.snaqr and self.snaqm:
			self.display_help("You can't use multiple presets!")
		elif self.snaq and self.snaqr:
			self.display_help("You can't use multiple presets!")
		elif self.snaq and self.snaqm:
			self.display_help("You can't use multiple presets!")
		elif self.snaqr and self.snaqm:
			self.display_help("You can't use multiple presets!")



		


		if not self.loci and not self.nex and not self.fas:
			print("Warning: No output type chosen. Are you sure you wanted this?? filterLoci.py will only report number of loci passing filters.")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nfilterLoci.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-l /path/to/loci\n")
		print ("Description: Filters pyRAD .loci file for N SNPs and N samples")

		print("\nProbably will add more functions later...")

		print("""
	Arguments:
		INPUT FILES
		-i,--input  : Input file as PHYLIP
		-P,--popmap : (Optional) Population map input file
		-w	: (Optional) List of white-listed loci (integers)
		PARAMETERS
		-s,--samples	: Minimum number of samples to retain locus [def=1]
		-p,--pis	: Minimum number of parsimony-informative sites [def=1]
		-n,--nex	: Boolean. Write loci as individual NEXUS files [false]
		-f,--fas	: Boolean. Write loci as individual FASTA files [false]
		-r,--rand   : Sample <x> random passing loci 
		-o,--out	: Output file prefix

		PRESETS
		-S,--snaq   : (only if -P) Preps .loci file for my SNaQ pipeline
			--Subsamples all --popmap populations to 1 sample per pop
			--Finds the individual with highest coverage and keeps it
			--Only retains loci containing 1 sample per popmap pop
		-R,--snaqR  : --snaq but randomly sampling one sample per pop
		-M,--snaqM  : --snaq but allows individuals to be missing from a locus
		-V,--pomo   : Subsamples UP TO [x] individuals per pop

""")
		sys.exit()


#Call main function
if __name__ == '__main__':
	main()
