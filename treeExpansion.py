#!/usr/bin/python

import re
import sys
import os
import getopt


def main():
	params = parseArgs()

	popsList = dict()
	#parse popmap file for dictionary of sample assignments
	if params.popmap:
		#print("Parsing popmap file...")
		popsList = parsePopmap_alt(params.popmap)

		if not params.tree:

			print("ERROR: No tree provided.")
			sys.exit(1)


		newtree = params.tree
		for pop in popsList:
			replace = ", ".join(popsList[pop])
			newtree = newtree.replace(str(pop), str(replace))
		print(newtree)


	else:
		print("ERROR: Popmap file must be provided.")
		sys.exit(1)


#function reads a tab-delimited popmap file and return dictionary of assignments
#function returns dict of pops, each pointint to list of taxa
def parsePopmap_alt(popmap):
	ret = dict()
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
					if len(stuff)!= 2:
						print("Uh oh! Record missing a field: ",stuff)
						continue
					if stuff[1] not in ret:
						l = list()
						l.append(stuff[0])
						ret[stuff[1]] = l
					else:
						ret[stuff[1]].append(stuff[0])
			return(ret)
		except IOError as e:
			print("Could not read file %s: %s"%(popmap,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(popmap,e))
			sys.exit(1)
		finally:
			fh.close()

#function returns first readable line from a file
#good for getting headers etc
def firstLine(f):
	with open(f, 'r') as fh:
		try:
			for line in fh:
				line = line.strip()
				if not line:
					continue
				else:
					return(line) #returns first real line
		except IOError as e:
			print("Could not read file %s: %s"%(f,e))
			sys.exit(1)
		except Exception as e:
			print("Unexpected error reading file %s: %s"%(f,e))
			sys.exit(1)
		finally:
			fh.close()


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 't:s:p:h', \
			["tree=","popmap="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.tree=None
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
			if opt in ('t','tree'):
				self.tree = firstLine(arg)
			elif opt in ('p', 'popmap'):
				self.popmap = arg
			elif opt in ('h', 'help'):
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.tree:
			self.display_help("Error: Missing required tree (--tree or --stree)")
		if not self.popmap:
			self.display_help("Error: Missing required popmap file (-p, --popmap)")


	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\ntreeExpansion.py\n")
		print ("Contact:Tyler K. Chafin, University of Arkansas,tkchafin@uark.edu")
		print ("Description: Expands Newick tree of clades to include all taxa in a popmap file")

		print("""
	Arguments:
		-p,--popmap	: Tab-delimited population map
		-t,--tree	: Newick tree in a file
		    or
		-s,--stree	: Newick tree given as a string
		-h,--help	: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
