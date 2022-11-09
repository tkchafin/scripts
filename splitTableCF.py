#!/usr/bin/python

import sys
import os
import warnings
import getopt
import toytree as tt
import pandas as pd
import numpy as np

warnings.simplefilter(action='ignore', category=FutureWarning)

def main():
	params = parseArgs()
	
	#read data
	with open(params.samples) as f:
		l = f.read().splitlines()
	cf = pd.read_csv(params.cf, header=0)
	tree = tt.tree(params.tree)
	
	#calculate mean ngenes for each sample
	cov=dict()
	for f in l:
		b=cf.eq(f).any(1)
		cov[f] = np.mean(cf[b]["ngenes"])
	
	#find which sample has best representation; will be used as placehold for whole list of samples
	placeholder = max(cov, key = cov.get)
	
	#make subset datasets
	removes = [s for s in l if s != placeholder]
	left_cf = subset_df_blacklist(cf, removes) #keeps placeholder
	removes2 = [s for s in tree.get_tip_labels() if s not in l]
	right_cf = subset_df_blacklist(cf, removes2)
	left_tree = tree_remove_blacklist(tree, removes)
	#print(left_tree.get_tip_labels())
	right_tree = tree_remove_whitelist(tree, l)
	#print(right_tree.get_tip_labels())
	
	#write outputs
	#ingroup_tree
	right_tree.write("ingroup_tree.tre", tree_format=5)
	#ingroup_cfs
	right_cf.to_csv("ingroup_cfs.csv", index=False, index_label=False)
	#outgroup_tree
	left_tree.write("outgroup_tree.tre", tree_format=5)
	#outgroup_cfs
	left_cf.to_csv("outgroup_cfs.csv", index=False, index_label=False)

def tree_remove_whitelist(tree, goodbois):
	all_tips = tree.get_tip_labels()
	rem = [a for a in all_tips if a not in goodbois]
	return(tree.drop_tips(names=rem))

def tree_remove_blacklist(tree, badbois):
	all_tips = tree.get_tip_labels()
	rem = [r for r in badbois if r in all_tips]
	return(tree.drop_tips(names=rem))

def subset_df_blacklist(df, badbois):
	ret = df.copy()
	for i in badbois:
		bools = ret.eq(i).any(1)
		ret = ret[~bools]
	return(ret)


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hc:s:n:t:', \
			["help", "cf=", "name=", "samples=", "tree="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.cf = None
		self.samples=None
		self.write="both"
		self.tree=None

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
			elif opt=="cf" or opt=="c":
				self.cf=arg
			elif opt=="samples" or opt=="s":
				self.samples=arg
			elif opt=="tree" or opt=="t":
				self.tree=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.tree or not self.samples or not self.cf:
			self.display_help("No files provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nsplitTableCF.py\n")
		print("Author: Tyler K Chafin, University of Colorado")
		print ("Contact: tyler.chafin@colorado.edu")
		print ("Description: Subsets a TableCF file (PhyloNetworks) given a list of samples comprising a monophyletic clade -- right now only designed for 1 split at a time")
		print("""
		-c,--cf		: CF table
		-s,--samples	: File with list of samples
		-t,--tree	: Tree file 
		-h,--help	: Help menu
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
