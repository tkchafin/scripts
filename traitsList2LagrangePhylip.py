#!/usr/bin/python

import sys
import os
import getopt

def main():
	params = parseArgs()
	
	traits=set()
	samples=dict()
	map=dict()
	
	if params.map:
		with open(params.map, "r") as m:
			for line in m:
				line=line.strip()
				if len(line)==0:
					continue
				stuff=line.split()
				if len(stuff) >2: 
					print("ERROR: Too many elements --",line)
				else:
					if stuff[0] in map.keys():
						map[stuff[0]].append(stuff[1])
					else:
						map[stuff[0]]=list()
						map[stuff[0]].append(stuff[1])

	with open(params.tab, "r") as t:
		for line in t:
			line=line.strip()
			if len(line) == 0:
				continue
			stuff=line.split("\t")
			if len(stuff) >2: 
				print("ERROR: Too many elements --",line)
			else:
				samples[stuff[0]]=set()
				#print(stuff[0])
				#print(stuff)
				if len(stuff) >1:
					splitstuff=stuff[1].split(",")
					for s in splitstuff:
						loc=s
						if params.map and s in map.keys():
							loc=map[s]
							for l in loc:
								samples[stuff[0]].add(l)
								traits.add(l)
							continue
						else:
							samples[stuff[0]].add(loc)
							traits.add(loc)
				#print(samples)
				#sys.exit()
		t.close()
	#print(traits)
	#sys.exit()


	trlen=len(traits)
	slen=len(samples)
	output=""
	rep=False
	for samp in samples:
		#print(samples[samp])
		oline = str(samp) + "\t"
		#if no traits, report and skip
		#print(samples[samp])
		if len(samples[samp]) < 1:
			if not rep:
				rep=True
				print("Samples were found without any trait data. Skipping samples:")
			print(samp)
			slen-=1
			continue
		else:
			for t in traits:
				if t in samples[samp]:
					oline = oline + "1"
				else:
					oline = oline + "0"
			oline+="\n"
			output = output+oline
		#sys.exit()
	#print(output)
	#write lagrange phylip file 
	with open(params.out, "w") as ofh:
		header=str(slen) + "\t" + str(trlen) + "\t(" + str(" ".join(traits)) + ")\n"
		print("Traits output in this order:")
		print(str(", ".join(traits)))
		ofh.write(header)
		ofh.write(output)
		ofh.close()
	

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'ht:o:m:', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.tab=None
		self.out="out.phy"
		self.map=None


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
			elif opt == "t":
				self.tab=arg
			elif opt=="o":
				self.out=arg
			elif opt=="m":
				self.map=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.tab:
			self.display_help("No table provided.")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\ntraitsList2LagrangePhylip.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Converts table of the form Sample \t Trait,Trait,Trait to phylip 0/1 format, for LAGRANGE of BioGeoBEARS")
		print("""
		-t:	Tab-delimited trait table
		-m: Option tab-delimited map grouping trait names
		-o: Output file name [default=out.phy]
""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
