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
	d=parseBarcodes(df)
	printBarcodes(d)

def parseBarcodes(df):

	set1 = df.iloc[:,0:6]
	set2 = df.iloc[:,7:13]
	set3 = df.iloc[:,14:20]

	d = dict()
	for part in [set1,set2,set3]:
		block = False
		count = 1
		for row_index,row in part.iterrows():
			#If block starts new specimen block
			if row[0] == "Specimen":
				block = True
				continue
			if row[0] in ["nan","NaN"]: #blank line, end of block
				block = False
				continue
			if block == True:
				if count <= 8:
					#print("sample:",row[0])
					d[row[0]]=row[-1]
					count += 1
				else:
					count = 1
					block = False
	return(d)

def printBarcodes(d):
	for key,val in d.items():
		row = key + "\t" + val
		print(row)


#Call main function
if __name__ == '__main__':
	try:
		main()
	except KeyboardInterrupt:
		sys.exit(1)
