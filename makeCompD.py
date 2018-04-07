#!/usr/bin/python

import sys
import os
import getopt
import itertools
import random

def main():
	params = parseArgs()

	if params.infile:
		pop_assign = dict()
		pops = dict() #Dict of lists containing samples per pop

		#parse popmap file for dictionary of sample assignments
		if params.popmap:
			print("\nParsing popmap file...")
			pop_assign = parsePopmap(params.popmap)
			pops = popmap2popdict(pop_assign)
			print("\tFound %s different populations and %s individuals!"%(len(pops), len(pop_assign)))
			#print(pops)
		else:
			print("ERROR: Popmap file must be provided.")
			sys.exit(1)

		print("\nParsing input tests file...")
		tests = parseDtests(params.infile)
		tl = len(tests)

		#Make sure proper number of lines
		if tl == 4:
			print("\t4 lines found in %s! Looks like you are doing the 4-taxon D-test?"%params.infile)
		elif tl == 5:
			print("\t5 lines found in %s! Looks like you are doing the Dfoil or partitioned-D?"%params.infile)
		else:
			print("\t%s lines found in %s! Is something wrong?"%(tl, params.infile))

		#Calulate number of tests
		num = combinations2Dlist(tests)
		if num in [None, 0]:
			print("\tUh oh! Calculated 0 tests. Something went wrong. Sorry I don't have better error checking in this script!")
		else:
			print("\tCalculated a total of %s separate D-tests! Will write %s files ending with %s!"%(num, num, params.out))

		#Get all tests
		split_tests = list(itertools.product(*tests))

		#If stats file provided, parse it
		weights = dict()
		if params.stats:
			print("\nFound a file containing individual weights! Parsing",params.stats)
			weights = parseWeights(params.stats)
			print("\tFound %s samples with weight data. Any samples missing data will not be selected."%len(weights))


		#For each test, write a Comp-D testfile
		print("\nBuilding individual Comp-D test files...")
		rand = True
		if params.stats:
			rand = False
			print("\t%s samples with highest weight from %s will be chosen per population"%(params.max, params.stats))
		else:
			print("\t%s samples per population will be selected at random."%params.max)
		print("\tWorking on %s tests..."%num)

		#For each test, build output file
		total_tests = 0
		passed_tests = 0
		for t in split_tests:
			outname = None
			chosen_inds= list() #list of lists will hold output samples
			#For each pop, grab individuals
			passed = True
			for pop in t:
				if outname:
					outname = outname + "+" + pop
				else:
					outname = pop

				if pop not in pops:
					passed = False

				if passed == False:
					continue
				#Select individuals for the pop
				inds = list()
				if rand:
					inds = chooseIndividualsRandom(pops[pop], params.max)
				else:
					inds = chooseIndividualsWeighted(pops[pop], weights, params.max)
					#print(inds)
				if not inds:
					passed = False
					continue
				if len(inds) <= 0:
					passed = False
					continue
				else:
					chosen_inds.append(inds)

			if passed == True:
				outname = outname + "." + params.out
				total_tests = total_tests + combinations2Dlist(chosen_inds)
				passed_tests = passed_tests + 1
				#Write output file
				with open(outname, 'w') as fh:
					try:
						for taxon in chosen_inds:
							l = ""
							for ind in taxon:
								l = l + " " + ind
							l = l + "\n"
							fh.write(l)
					except IOError:
						print("Could not open file ",outname)
						sys.exit(1)
					finally:
						fh.close()
			else:
				print("\tUh oh! Test %s failed! Skipping it."%outname)

		print("\tSuccessfully wrote %s test files, with an average of %s replicates per test!"%(passed_tests, total_tests/passed_tests))
		print("\tDone!\n")

#Function to return random list of max_num individuals from a list
def chooseIndividualsRandom(l, max_num):
	max_num = int(max_num)
	if len(l) <= max_num:
		return(l)
	else:
		return(random.sample(l, max_num))

#Function to return random list of max_num individuals from a list
def chooseIndividualsWeighted(l, w, max_num):
	max_num = int(max_num)
	if len(l) <= max_num:
		return(l)
	else:
		subset = dict()
		for ind in l:
			if ind in w:
				subset[ind] = w[ind]
			else:
				subset[ind] = 0
		_sorted = sorted(subset, key=subset.get, reverse=True)
		return(_sorted[0:max_num])


#Function to calculate number of combinations from a 2D list
def combinations2Dlist(ll):
	total = None
	for l in ll:
		if total:
			total = total * len(l)
		else:
			total = len(l)
	return(total)

#Function to grab weights and return dict of individuals with weight
def parseWeights(w):
	if os.path.exists(w):
		weights = dict()
		with open(w, 'r') as fh:
			try:
				lines = []
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					l = line.split()
					weights[l[0]] = int(l[1])
				#print(weights); ;sys.exit()
				return(weights)
			except IOError:
				print("Could not read file ",w)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%w)


#Parse D tests file
def parseDtests(t):
	if os.path.exists(t):
		with open(t, 'r') as fh:
			try:
				lines = []
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					l = line.split()
					lines.append(l)
				return(lines)
			except IOError:
				print("Could not read file ",t)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%t)

#Function to parse a popmap (dict) and return dict of lists (list per pop)
def popmap2popdict(popmap):
	p = dict()
	for ind in popmap:
		if popmap[ind] in p:
			p[popmap[ind]].append(ind)
		else:
			l = [ind]
			p[popmap[ind]] = l
	return(p)


#Function to check that list of sample names and popmap entries match
def validatePopmap(samples, popmap):
	print(samples)
	print(popmap)
	for samp in samples:
		if samp in popmap:
			print("Warning: Sample %s not found in popmap!"%samp)
	for key in popmap:
		if key not in samples:
			print("Warning: Sample %s found in popmap has no data!"%key)

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


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'p:w:hi:o:m:', \
			["popmap=","help","weights=","in=","out=","max="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.popmap=None
		self.out="compd"
		self.stats=None
		self.infile=None
		self.max=10

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
			if opt in ('p', 'popmap'):
				self.popmap = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('o','out'):
				self.out = arg
			elif opt in ("w","weights"):
				self.stats = arg
			elif opt in ("i","in"):
				self.infile = arg
			elif opt in ("m","max"):
				self.max = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.popmap:
			self.display_help("Error: Missing required popmap file (-p,--popmap).")
		if not self.infile:
			self.display_help("Error: Missing required test file (-i,--in).")



	def display_help(self, message=None):
		if message is not None:
			print("")
			print (message)
		print ("\nmakeCompD.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-p </path/to/popmap> -i </path/to/infile> <-s /path/to/stats> <-m int>\n")
		print ("Description: Creates input test files for Comp-D (https://github.com/smussmann82/Comp-D_MPI)")

		print("""
	Arguments:
		-p,--popmap	: Path to tab-delimited population map
			-Assumes that all sampled in popmap can be used for CompD tests
		-i,--infile	: Prefix for input file
			Using format of CompD, but with Population IDs from the popmap file
		-w,--weights	: Optional. File containing individual weights
			-Should be formatted as sampleName \\t weight
			-Samples with highest weight are chosen. For example, this could be
			 provided as number of loci present. Samples with most loci will be selected
			-Samples missing from weights file will be assigned weights of 0
		-m,--max	: Maximum individuals to sample per populations <default=10>
			-If used without <-s>, samples will be randomly chosen
			-If used with <-s>, x samples with highest weight are chosen
		-o,--out	: Suffix for output file
			-Prefix will be formatted as: O+P3+P2+P1.$out for 4-taxon test
			-$out be default is 'compd'
		-h,--help	: Displays help menu
""")
		print("""
	Infile format:
		OUT1 OUT2 OUT3+OUT4
		P3A P3B
		P2A P2B P2C P2D P2D
		P1A P1B

	#Order follows that of Comp-D
	#Pop IDs will be used to create CompD infiles for all combinations

	""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
