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
		pops = dict()
		enumCounter = 1
		gen_cats = 0

		numPops = 0

		print("\nGetting assignments of individuals...")
		for key, value in probs.items():
			index += 1;
			if ("nan" in value):
				nan += 1;
				continue #skip individuals which couldn't be assigned

			if p[index] not in pops:
				#add pop to pops
				numPops += 1
				pops[p[index]] = [0,0,0,0,0,0,0] #P0,P1,F1,F2,B0,B1,FN

			assign=False
			for i,v in enumerate(value):
				if float(v) >= params.thresh:
					assign=True
					pops[p[index]][i] += 1
			if assign == False:
				pops[p[index]][6] += 1

		#Convert pops counts to proportions
		pop_props = dict()

		print("\nWriting INDIVQ file for distruct (with pops as individuals): PROPS_indivq.txt")

		with open ("PROPS_indivq.txt", "w") as IQ:
			try:
				popNum=0
				for key, value in pops.items():
					if key not in pop_props:
						pop_props[key] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
					total = 0
					for count in value:
						total += count

					for i,v in enumerate(value):
						pop_props[key][i] = str("{0:.4f}".format((v/total)))

					#print(key,":","\t".join(pop_props[key]))
					indline = str(popNum) + "\t" + str(popNum) + "\t(0)\t" + str(popNum+1) + "\t: " + "\t".join(pop_props[key]) + "\n"
					IQ.write(indline)
					popNum+=1
			except IOError as e:
				print("Could not open file:",e)
				sys.exit(1)
			finally:
				IQ.close()

		print("Writing PROPS TABLE file for reading in R : PROPS_table.tsv")

		with open ("PROPS_table.tsv", "w") as IQ:
			try:
				indline = "POP\tP0\tP1\tF1\tF2\tB0\tB1\tFN\n"
				IQ.write(indline)
				for key, value in pops.items():
					if key not in pop_props:
						pop_props[key] = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
					total = 0
					for count in value:
						total += count

					for i,v in enumerate(value):
						pop_props[key][i] = str("{0:.4f}".format((v/total)))

					#print(key,":","\t".join(pop_props[key]))
					indline = str(key) + "\t" + "\t".join(pop_props[key]) + "\n"
					IQ.write(indline)
			except IOError as e:
				print("Could not open file:",e)
				sys.exit(1)
			finally:
				IQ.close()

		print("Writing dummy POPQ file: PROPS_popq.txt")
		with open("PROPS_popq.txt", "w") as PQ:
			try:
				for i in range(0,numPops):
					numPop = i + 1
					line = str(numPop) + ":\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t1\n"
					PQ.write(line)
			except IOError:
				print("Could not open file",params.popq)
				sys.exit(1)
			finally:
				PQ.close()

		print("Writing Labels file: PROPS_labels.txt")
		with open("PROPS_labels.txt", "w") as ID:
			try:
				popNum = 0
				for pop, enum in pops.items():
					popNum += 1
					out = str(popNum) + " " + str(pop) + "\n"
					ID.write(out)

			except IOError:
				print("Could not open file PROPS_labels.txt")
				sys.exit(1)
			finally:
				ID.close()

		print("Writing COLOR PERMUTATION file: PROPS_geno.perm")
		with open("PROPS_geno.perm", "w") as PERM:
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

				FN = "7 brown\n"
				print("FN: Brown (brown)")
				PERM.write(FN)

				print()

			except IOError:
				print("Could not open file PROPS_geno.perm")
				sys.exit(1)
			finally:
				PERM.close()

		print("Writing Distruct paramsfile: PROPS_params.txt")
		with open("PROPS_params.txt", "w") as PAR:
			try:
				stuff = getParams(numPops)
				PAR.write(stuff)

			except IOError:
				print("Could not open file PROPS_params.txt")
				sys.exit(1)
			finally:
				PAR.close()


		print("\nDone!\n")
	else:
		print("Missing required inputs.")
		sys.exit(1)

def getParams(np):
	par = """
#define INFILE_POPQ PROPS_popq.txt
#define INFILE_INDIVQ PROPS_indivq.txt
#define INFILE_LABEL_BELOW PROPS_labels.txt
#define INFILE_LABEL_ATOP PROPS_labels.txt
#define INFILE_CLUST_PERM PROPS_geno.perm
#define OUTFILE PROPS.ps
#define K 7
"""
	par = par + "#define NUMPOPS " + str(np) + "\n"
	par = par + "#define NUMINDS " + str(np) + "\n"
	par = par + """#define PRINT_INDIVS 1
#define PRINT_LABEL_ATOP 1
#define PRINT_LABEL_BELOW 0
#define PRINT_SEP 1
#define FONTHEIGHT 6
#define DIST_ABOVE -160
#define DIST_BELOW -50
#define BOXHEIGHT 150
#define INDIVWIDTH 20
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
			options, remainder = getopt.getopt(sys.argv[1:], 'i:p:t:h', \
			["pops=","input=",'thresh=',"help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.pops = None
		self.pofz=None
		self.out = "NH_indivq.txt"
		self.popq = "NH_popq.txt"
		self.thresh=float(0.80)

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
			elif opt in ('t', 'thresh'):
				self.thresh = float(arg)
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.pops:
			self.display_help("Error: Missing required PopID file (-p, --pops)")
		if not self.pofz:
			self.display_help("Error: Missing required PofZ file (-i, --input)")

		if self.thresh <= 0.50:
			self.display_help("<-t,--thresh> cannot be less than 0.50")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nnewhybs2props.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("\nUsage: ", sys.argv[0], "-i aa-PofZ.txt -p popmap \n")
		print ("Description: Creates proportion tables which can be plotted in DISTRUCT or R")

		print("""
	Arguments:
		-i,--input	: aa-PofZ.txt output from NewHybrids.
		-p,--pops	: Path to population IDs for NewHybrids samples
			Format: List of population IDs in the SAME ORDER as NewHybrids output.
			Note: My phylip2newhybrids.pl script will create this for you.
		-t,--thresh	: Threshold posterior probability to assign sample to class [default=0.8]
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
