#!/usr/bin/python

import sys
import pandas as pd

def main():
	x=""
	if sys.argv[1]:
		#print(sys.argv[1])
		x=sys.argv[1]
	else:
		print("Exiting because .xlsl not provided")
		print("Usage:",sys.argv[0]," <L####-##.xlsx>")
		sys.exit(1)

	df = pd.read_excel(x,sheet_name=0)
	#print(df)
	parseDnaex(df)


def parseDnaex(df):
	block = False
	new = pd.DataFrame(columns=["sample","field","concentration"])
	for index,row in df.iterrows():
		#If block starts new specimen block
		if str(row[0]) == "spec.":
			block = True
			continue
		if str(row[0]) in ["nan","NaN"]: #blank line, end of block
			continue
		if block == True:
			if str(row[1]) in ["nan","NaN"]: #blank line, end of block
				continue
			name = str(row[1]) + str(row[2]) + str(row[3])
			field = str(row[4])
			conc = str(row[-1])
			print(name, end="\t")
			print(field,end="\t")
			print(conc,end="\n")


#Call main function
if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.exit(1)
