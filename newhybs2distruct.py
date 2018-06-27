#!/usr/bin/python

import sys
import os
import getopt
import collections

def main():
	params = parseArgs()
	if params.pops and params.pofz:
		#get pop IDS
		p = readList(params.pops)

		#get OrderedDict of prob results
		probs = readNewHybs(params.pofz)

		#iterate over probs to make output
		nan = 0
		index = -1
		popCount = dict()
		popEnum = dict()
		enumCounter = 1
		gen_cats = 0

		numPops = 0
		numInds = 0
		print("\nWriting INDIVQ file for distruct:",params.out)

		with open (params.out, "w") as IQ:
			try:
				for key, value in probs.items():
					index += 1;
					if ("nan" in value):
						nan += 1;
						continue #skip individuals which couldn't be assigned

					numInds += 1
					#track population ID and count per pop
					if (p[index] not in popEnum):
						popEnum[p[index]] = enumCounter
						enumCounter += 1
						numPops += 1
					if (p[index] not in popCount):
						popCount[p[index]] = 1
					else:
						popCount[p[index]] += 1

					if gen_cats == 0:
						gen_cats = len(value)
					elif gen_cats != len(value):
						print("Warning: Samples don't have the same number of probabilities! Something is wrong")

					#build output line for INDIVQ
					indline = str(key) + "\t" + str(key) + "\t(0)\t" + str(popEnum[p[index]]) + "\t: " + "\t".join(value) + "\n"
					IQ.write(indline)
					#print(key, "(",p[index], "): ", value)

				if (nan > 0):
					print("Warning:",nan,"individuals had \"nan\" probabilities are were skipped.")
			except IOError:
				print("Could not open file",params.out)
				sys.exit(1)
			finally:
				IQ.close()

		print("Writing dummy POPQ file:", params.popq)
		with open(params.popq, "w") as PQ:
			try:
				for pop, enum in popEnum.items():
					out = str(enum) + ":"
					for cat in range(gen_cats):
						out = out + "\t0.0"
					out = out + "\t" + str(popCount[pop]) + "\n"
					PQ.write(out)

			except IOError:
				print("Could not open file",params.popq)
				sys.exit(1)
			finally:
				PQ.close()

		print("Writing Labels file: NH_labels.txt")
		with open("NH_labels.txt", "w") as ID:
			try:
				for pop, enum in popEnum.items():
					out = str(enum) + " " + str(pop) + "\n"
					ID.write(out)

			except IOError:
				print("Could not open file NH_labels.txt")
				sys.exit(1)
			finally:
				ID.close()

		print("Writing COLOR PERMUTATION file: NH_geno.perm")
		with open("NH_geno.perm", "w") as PERM:
			try:
				print()
				P1 = "1 RdGy_6_div_1\n"
				print("P1: Red (RdGy_6_div_1)")
				PERM.write(P1)

				P2 = "2 RdBu_6_div_6\n"
				print("P2: Blue (RdBu_6_div_6)")
				PERM.write(P2)

				F1 = "3 Greens_6_seq_5\n"
				print("F1: Green (Greens_6_seq_5)")
				PERM.write(F1)

				F2 = "4 Greens_6_seq_2\n"
				print("F2: Light Green (Greens_6_seq_2)")
				PERM.write(F2)

				BO1 = "5 RdBu_6_div_3\n"
				print("BO1: Light red (RdBu_6_div_3)")
				PERM.write(BO1)

				BO2 = "6 RdBu_6_div_4\n"
				print("BO2: Light Blue (RdBu_6_div_4)")
				PERM.write(BO2)

				print()

			except IOError:
				print("Could not open file NH_geno.perm")
				sys.exit(1)
			finally:
				PERM.close()

		print("Writing Distruct paramsfile: NH_params.txt")
		with open("NH_params.txt", "w") as PAR:
			try:
				stuff = getParams(numPops, numInds)
				PAR.write(stuff)

			except IOError:
				print("Could not open file NH_params.txt")
				sys.exit(1)
			finally:
				PAR.close()


		print("Done!\n")
	else:
		print("Missing required inputs.")
		sys.exit(1)

def getParams(np, ni):
	par = """
#define INFILE_POPQ NH_popq.txt
#define INFILE_INDIVQ NH_indivq.txt
#define INFILE_LABEL_BELOW NH_labels.txt
#define INFILE_LABEL_ATOP NH_labels.txt
#define INFILE_CLUST_PERM NH_geno.perm
#define OUTFILE NH.ps
#define K 6
"""
	par = par + "#define NUMPOPS " + str(np) + "\n"
	par = par + "#define NUMINDS " + str(ni) + "\n"
	par = par + """#define PRINT_INDIVS 1
#define PRINT_LABEL_ATOP 1
#define PRINT_LABEL_BELOW 0
#define PRINT_SEP 1
#define FONTHEIGHT 6
#define DIST_ABOVE -160
#define DIST_BELOW -50
#define BOXHEIGHT 150
#define INDIVWIDTH 2
#define ORIENTATION 1
#define XORIGIN 200
#define YORIGIN 10
#define XSCALE 1
#define YSCALE 1
#define ANGLE_LABEL_ATOP 270
#define ANGLE_LABEL_BELOW 270
#define LINEWIDTH_RIM 3
#define LINEWIDTH_SEP 1
#define LINEWIDTH_IND 3
#define GRAYSCALE 0
#define ECHO_DATA 1
#define REPRINT_DATA 1
#define PRINT_INFILE_NAME 0
#define PRINT_COLOR_BREWER 1"""
	return(par)

#reads and returns a list from a file
def readList(l):
	if os.path.exists(l):
		with open(l, 'r') as fh:
			try:
				ret = list()
				for line in fh:
					line = line.strip()
					if not line:
						continue
					ret.append(line)
				return(ret)
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)


#reads assignment probabilities from NewHybs PofZ output file
def readNewHybs(p):
	if os.path.exists(p):
		with open(p, 'r') as fh:
			try:
				ret = collections.OrderedDict()
				count = 0;
				for line in fh:
					line = line.strip()
					if not line:
						continue
					count += 1
					if count == 1:
						continue #skip first non-blank line, which is the header
					else:
						arr = line.split()
						ret[arr[0]] = list(arr[2:])
				return(ret)
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
			options, remainder = getopt.getopt(sys.argv[1:], 'i:p:', \
			["pops=","input="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.pops = None
		self.pofz=None
		self.out = "NH_indivq.txt"
		self.popq = "NH_popq.txt"

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
			if opt in ('p', 'pops'):
				self.pops = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('i', 'input'):
				self.pofz = arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.pops:
			self.display_help("Error: Missing required PopID file (-p, --pops)")
		if not self.pofz:
			self.display_help("Error: Missing required PofZ file (-i, --input)")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nnewhybs2distruct.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-i aa-PofZ.txt -p popmap \n")
		print ("Description: Creates inputs for DISTRUCT from NewHybrids output.")

		print("""
	Arguments:
		-i,--input	: aa-PofZ.txt output from NewHybrids.
		-p,--pops	: Path to population IDs for NewHybrids samples
			Format: List of population IDs in the SAME ORDER as NewHybrids output.
			Note: My phylip2newhybrids.pl script will create this for you.
		-o,--out	: Output file name <default = out.nex>
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
