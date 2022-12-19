#!/usr/bin/python

import sys
import os
import getopt
#import toytree as tt
import random

def main():
	params = parseArgs()

	seqs=dict()	
	for f in read_phylip(params.phylip):
                        seqs[f[0]] = f[1]

	#tree=tt.tree(params.tree, tree_format=0)

	if not params.samples:
		params.samples=int(params.freq*len(list(seqs.keys())))
	
	print("Generating",params.reps,"random subsets of",params.samples,"samples eac")

	for r in range(params.reps):
		print("starting replicate",str(r))
		prefix=params.out + "_" + str(r)
		print("subsetting alignment")
		keeps=dict(random.sample(seqs.items(), params.samples))
		bad_bois=[k for k in seqs.keys() if k not in keeps]
		#print("subsetting tree")
		#stree = tree.drop_tips(names=bad_bois)
		print("writing subset file")
		write_phylip(prefix+".phylip",keeps)
		#stree.write(prefix+".tre", tree_format=0)


#Print dict to phylip file
def write_phylip(p, aln):
        with open(p, 'w') as fh:
                try:
                        header = getPhylipHeader(aln) + "\n"
                        fh.write(header)

                        for sample in aln.keys():
                                line = str(sample) + "\t" + "".join(aln[sample]) + "\n"
                                fh.write(line)
                except IOError as e:
                        print("Could not read file %s: %s"%(p,e))
                        sys.exit(1)
                except Exception as e:
                        print("Unexpected error reading file %s: %s"%(p,e))
                        sys.exit(1)
                finally:
                        fh.close()	

#Returns header for Phylip file from a dictionary of samples w/ data
def getPhylipHeader(d):
        numSamp = 0
        numLoci = None
        for sample in d:
                numSamp = numSamp + 1
                if not numLoci:
                        numLoci = len(d[sample])
                else:
                        if numLoci != len(d[sample]):
                                print("getPhylipHeader: Warning: Sequences of unequal length.")
        header = str(numSamp) + " " + str(numLoci)
        if numLoci == 0 or not numLoci:
                print("getPhylipHeader: Warning: No loci in dictionary.")
        if numSamp == 0:
                print("getPhylipHeader: Warning: No samples in dictionary.")
        return(header)


#Read samples as PHYLIP. Generator function
def read_phylip(phy):
        if os.path.exists(phy):
                with open(phy, 'r') as fh:
                        try:
                                num=0
                                for line in fh:
                                        line = line.strip()
                                        if not line:
                                                continue
                                        num += 1
                                        if num == 1:
                                                continue
                                        arr = line.split()
                                        yield(arr[0], arr[1])
                        except IOError:
                                print("Could not read file ",phy)
                                sys.exit(1)
                        finally:
                                fh.close()
        else:
                raise FileNotFoundError("File %s not found!"%phy)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hs:f:r:p:o:m:', \
			["help", "reps=","phylip=","out=", "method=", "samples=", "freq="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		#self.tree=None
		self.reps=10
		self.freq=0.1
		self.samples=None
		self.phylip=None
		self.method="random"
		self.out="subset"


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
			elif opt=="phylip" or opt=="p":
				self.phylip=arg
			elif opt=="method" or opt=="m":
				self.method=arg
			elif opt=="reps" or opt=="r":
				self.reps=int(arg)
			elif opt=="freq" or opt=="f":
				self.freq=float(arg)
			elif opt=="samples" or opt=="s":
				self.samples=int(arg)
			elif opt=="out" or opt=="o":
				self.out=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.phylip and not self.tree:
			self.display_help("Must provide input tree (newick) and alignment (phylip) files.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nalignment_subsetter.py\n")
		print ("Description: Generate random subsets of an input phylip (alignment)")
		print("""
		-p,--phylip	: Path to input phylip file
		-s,--samples	: Number of samples to keep
		-f,--freq	: Sampling frequency (must set either -f or -s)
		-r,--reps	: Number of replicates to generate
		-o,--out	: Output file name (default=out.fas)
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
