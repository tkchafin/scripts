#!/usr/bin/python

import os
import sys 

def main():
	if len(sys.argv) <= 1:
		print("No files provided!")
		print('Usage: ./concatFasta.py *.fasta - or - ./concatFasta.py 1.fas 2.fas...')
		sys.exit(1)
	files = sys.argv[1:]

	print("concatenating fastas using the following order:")
	print("if this is incorrect, change something")

	pre=None
	samps=dict()
	for file in sorted(files):
		print(file)
		pre=file.split("_")[0]
		for s in read_fasta(file):
			if s[0] not in samps:
				samps[s[0]] = s[1]
			else:
				samps[s[0]] = samps[s[0]] + s[1]

	print("Using prefix from files to write output:",pre)
	oname = pre + ".fasta"
	write_fasta(oname, samps)

#write fasta from dict
def write_fasta(name, d):
	with open(name, 'w') as fh:
		try:
			for sample in d.keys():
				to_write = ">" + str(sample) + "\n" + d[sample] + "\n"
				fh.write(to_write)
		except IOError as e:
			print("Could not read file:",e)
			sys.exit(1)
		except Exception as e:
			print("Unexpected error:",e)
			sys.exit(1)
		finally:
			fh.close()

#Read samples as FASTA. Generator function
def read_fasta(fas):
	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#print(line)
					if line[0] == ">": #Found a header line
						#If we already loaded a contig, yield that contig and
						#start loading a new one
						if contig:
							yield([contig,seq]) #yield
							contig = "" #reset contig and seq
							seq = ""
						split_line = line.split()
						contig = (split_line[0].replace(">",""))
					else:
						seq += line
				#Iyield last sequence, if it has both a header and sequence
				if contig and seq:
					yield([contig,seq])
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)
		
		

#Call main function
if __name__ == '__main__':
    main()