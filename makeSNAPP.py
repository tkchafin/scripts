#!/usr/bin/python

import getopt
import sys
import os

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'x:p:l:o:h', \
			["xml=","phy=","list=","out=","help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.xml=None
		self.phy=None
		self.tax=None
		self.out="out.xml"

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			print(opt,arg)
			if opt in ('x', 'xml'):
				self.xml = arg
			elif opt in ('h', 'help'):
				pass
			elif opt in ('p','phy'):
				self.phy = arg
			elif opt in ('l','list'):
				self.tax = arg
			elif opt in ('o','out'):
				self.out = out
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		self.xml or self.display_help("INPUT ERROR: No XML provided")
		self.phy or self.display_help("INPUT ERROR: No PHYLIP provided")
		self.tax or self.display_help("INPUT ERROR: No TAXON LIST provided")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nmakeSNAPP.py\n")
		print ("Contact:\n\n\tTyler K. Chafin\n\tUniversity of Arkansas\n\ttkchafin@uark.edu\n")
		print ("\nUsage:\n\t", sys.argv[0], "-x </path/to/xml>\n")
		print ("Description:\n")
		print("\tmakeSNAPP.py is a quickly written and shitty script to help make lots\n",\
		"\tof SNAPP XML input files from a given list of taxa to include, a phylip file\n",\
		"\twith sequence data to pull from, and an axample XML file containing parameter \n",\
		"\tsetting and any priors. It was made for a very specific purpose. \n")

		print("""
		Input options:

			-x,--xml	: XML file created in Beauti
			-p,--phy	: Phylip file
			-l,--list	: .txt file containing a list of taxa to subset
			-o,--out	: (Optional) output prefix [default:out.xml]
			-h,--help	: Displays help menu""")
		print()
		sys.exit()


################################# MAIN #########################################
params = parseArgs()
