#!/usr/bin/python

import os
import sys
import argparse

def main():
		params = parseArgs()
		output=list()
		with open(params.input, "r") as fh:
			coord=1
			linecount=0
			chrom=None
			start=1
			step=1
			current_start=None
			report=True
			total=0
			for line in fh:
				line = line.strip()
				if not line:
					continue
				linecount += 1
				if linecount == 1 or "=" in line:
					header = line.split()
					for field in header:
						parts=field.split("=")
						if parts[0] == "chrom":
							if chrom and parts[1] != chrom:
								print("Found new chrom:",parts[1])
							chrom=parts[1]
						elif parts[0] == "start":
							start = int(parts[1])
						elif parts[0] == "step":
							step = int(parts[1])
						else:
							continue
					if current_start and start > coord:
						print("Started new region:",start, "- jumped from",coord)
						padded_start = current_start - params.padding
						if padded_start <= 0:
							padded_start = 1
						end = coord
						padded_end = end + params.padding
						if end - current_start > params.min_length:
							total += (padded_end-padded_start)
							oline=str(chrom)+":"+str(padded_start)+"-"+str(end)+"\n"
							output.append(oline)
						#print(oline)
						current_start = None
						coord=start
					else:
						coord=coord+step
					if not chrom:
						sys.exit("No chrom field found in header! Exiting script.")
					if report:
						print("\nChrom is:",chrom)
						print("Starting coordinate:",start)
						print("Step size:",step)
						print("Minimum phaseCons score:",params.min_score)
						print("Minimum length to report interval:",params.min_length)
						if params.padding > 0:
							print("Padding (+/-) for interval coordinates:",params.padding)
						print("\n--\n")
						report=False
					continue

				if float(line) >= params.min_score:
					if not current_start:
						current_start=coord
				else:
					#NOT above threshold. If there is a previous interval, check it now
					if current_start:
						padded_start = current_start - params.padding
						if padded_start <= 0:
							padded_start = 1
						end = coord
						padded_end = end + params.padding
						if end - current_start > params.min_length:
							#print(end-current_start+(2*params.padding))
							total += (padded_end-padded_start)
							oline=str(chrom)+":"+str(padded_start)+"-"+str(end)+"\n"
							output.append(oline)
							#print(oline)
						current_start = None
				coord = coord + step
		fh.close()

		print("\n--\nDone! Writing output to:", params.output)

		with open(params.output, "w") as ofh:
			if len(output) > 0:
				for l in output:
					ofh.write(l)
		ofh.close()
		print("\nProcess complete. Total bases included in retained intervals:",total, "\n")








#argument parsing
def parseArgs():
	help = """
	parsePhaseCons.py

	Author: Tyler K. Chafin
	Contact: tkchafin@uark.edu

	Description: Processes phaseCons outputs to generate a set of intervals with phaseCons score above X

	Input should be a file of phaseCons scores, with a header including the following information:
	chrom=<chrom_name> start=<start_coordinate> step=<step_size>
	"""
	parser = argparse.ArgumentParser(description=help)

	parser.add_argument('--min_length', dest='min_length', type=int, default=10,
						help='Minimum interval length to report [default=10]')
	parser.add_argument('--padding', dest='padding', type=int, default=0,
						help='Distance to pad interval coordinated (e.g. output=start-padding:end+padding) [default=0]')
	parser.add_argument('--min_score', dest='min_score', type=float, default=0.5,
						help='Minimum phaseCons score [default=0.5]')
	parser.add_argument('--input', dest='input', type=str,
						help='Input .pp.data file')
	parser.add_argument('--output', dest='output', type=str, default="phaseCons_intervals.bed",
						help='Output .bed file [default=phaseCons_intervals.bed]')

	args = parser.parse_args()

	if not args.input:
		sys.exit("Missing inputs")

	return args

#Call main function
if __name__ == '__main__':
    main()
