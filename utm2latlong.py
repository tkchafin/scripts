#!/usr/bin/python

import sys
import os
import getopt
import utm

def main():
	params = parseArgs()

	if params.utm:
		for line in readTSV(params.utm):
			if params.zone and params.hemi:
				coords = utm.to_latlon(float(line[1]), float(line[2]), params.zone, params.hemi)
				oline=str(line[0])+"\t"+str(coords[0])+"\t"+str(coords[1])
				print(oline)
			elif params.inline:
				z=line[3][:-1]
				h=line[3][-1]
				coords = utm.to_latlon(float(line[1]), float(line[2]), int(z), h)
				oline=str(line[0])+"\t"+str(coords[0])+"\t"+str(coords[1])
				print(oline)
			else:
				params.display_help("No UTM zone information provided.")
	elif params.latlong:
		for line in readTSV(params.latlong):
			coords = utm.from_latlon(float(line[1]), float(line[2]))
			oline=str(line[0])+"\t"+str(coords[0])+"\t"+str(coords[1]) + "\t" +str(coords[2])+str(coords[3])
			print(oline)
	else:
		params.display_help("No input provided")


#generator function, reads tsv line by line
def readTSV(tab):
	with open(tab, 'r') as fh:
		try:
			for line in fh:
				line = line.strip()
				if not line:
					continue
				yield(line.split())
		except IOError:
			print("Could not read file ",tab)
			sys.exit(1)
		finally:
			fh.close()

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'hu:c:z:l:i', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.utm=None
		self.latlong=None
		self.zone=None
		self.hemi=None
		self.inline=False


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
			elif opt == "c":
				self.latlong=arg
			elif opt=="l":
				self.hemi=arg
			elif opt == "z":
				self.zone=int(arg)
			elif opt=="u":
				self.utm=arg
			elif opt=="i":
				self.inline=True
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.utm and not self.latlong:
			self.display_help("No input file provided (must be one of: <-u> or <-c>)")
		if self.utm and self.latlong:
			self.display_help("Options not compatible: <-u> <-c>")
		if self.utm:
			if not self.zone and not self.hemi and not self.inline:
				self.display_help("Must provide zone number <-z> and letter <-l> with UTMs or as inline <-i>")
			if self.zone and self.inline:
				self.display_help("Options not compatible: <-i> <-z>")



	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\n<template.py>\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: ")
		print("""
		Arguments
		-u	: Tab-delimited table of UTM coordinates (2nd col = Easting; 3rd col = Northing)
		  -or-
		-c	: Tab-delimited table of lat/long coordinates (2nd col = lat; 3rd col= long)
		-z	: If converting UTM to lat/long, provide zone number here (e.g. "15")
		-l	: If converting UTM to lat/long, provide zone letter here (e.g. "N")
		-i	: If converting from UTMs, zone can be as 4th column (e.g. "15S")

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
