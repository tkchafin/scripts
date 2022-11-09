#!/usr/bin/python

import re
import sys
import os
import getopt
import sys
import pandas as pd
import io
import time
import pandas as pd
from io import StringIO
from Bio import Entrez
from Bio import AlignIO
from Bio import SeqIO
from Bio import Alphabet
from sqlite3 import OperationalError
from Bio.Align.Applications import MuscleCommandline


def main():
	params = parseArgs()

	#Specify email for Entrex session
	Entrez.email = params.email


	#get gi list for the query
	hits=list()
	if params.query:
		print( "\nQuerying NCBI GenBank database %s for query \"%s\"...\n" %(params.db, params.query))
		handle = Entrez.esearch( db=params.db,term=params.query,retmax=params.retmax )
		hits = Entrez.read(handle)['IdList']
	elif params.qlist:
		with open(params.qlist, 'r') as fh:
			try:
				for line in fh:
					line = line.strip()
					if not line:
						continue
					else:
						print( "\nQuerying NCBI GenBank database %s for query \"%s\"..." %(params.db, line))
						handle = Entrez.esearch( db=params.db,term=line,retmax=params.retmax )
						local_hits = Entrez.read(handle)['IdList']
						print("Found:",len(local_hits), "hits.")
						hits += local_hits
			except IOError as e:
				print("Could not read file:",e)
				sys.exit(1)
			except Exception as e:
				print("Unexpected error:",e)
				sys.exit(1)
			finally:
				fh.close()
	else:
		print("Error: No queries given.")
		sys.exit(1)

	#print(hits)

	#Report how many hits were found
	#probably could parse the hits and only DL ones we want,
	#but I don't care enough right now sorry
	if (len(hits) <= 0):
		print("No hits found! Try a better query. \n")
		sys.exit(0)
	print("Found %s entries! Downloading in batches of %s...\n"%(len(hits), params.batchSize))

	#If outfile desired, set up filehandle
	if params.Fonly:
		fetchToFile(hits, params.db, params.batchSize, params.out, "fasta")
	elif params.GBonly:
		fetchToFile(hits, params.db, params.batchSize, params.out, "gb")
	else:
		#post NCBI query
		search_handle     = Entrez.epost(params.db, id=",".join(hits))
		search_results    = Entrez.read(search_handle)
		webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]

		#fetch batches of results, parse each batch
		count=1

		#initialize table
		foundGenes = dict()
		coordinates = dict()
		accessions = dict()
		bad_accessions = dict()

		for start in range( 0,len(hits),params.batchSize ):

			print("Parsing batch %s..." %count)
			try:
				my_args = {'db':params.db, 'retmode':'xml', 'rettype' : params.rettype, 'retstart':start, 'retmax':params.batchSize, 'webenv':webenv, 'query_key':query_key}
				if params.api:
					my_args['api']=params.api
				handle = Entrez.efetch(**my_args)
				records = Entrez.parse(handle)

			except Exception as e:
				print("Oh no! Unknown exception on batch",count,":",e)
				count += 1
			else:
				#If no exception, parse records
				for rec in records:
					#print(rec)
					if (rec["GBSeq_moltype"] == "DNA"):
						try:
							gb_quals = rec["GBSeq_feature-table"][0]["GBFeature_quals"]
							gb_quals.extend(rec["GBSeq_feature-table"][1]["GBFeature_quals"])
							for qual in gb_quals:
								#print(rec)
								#print(qual)
								if qual["GBQualifier_name"] == "gene":
									#grab some information
									taxon = rec["GBSeq_organism"].strip().lower().replace(" ","_")
									gene = qual["GBQualifier_value"].strip().lower().replace(" ","_")
									seq = rec["GBSeq_sequence"].strip().lower()
									acc = rec["GBSeq_primary-accession"].strip()
									#print("Found %s (%s): %s\n"%(taxon, gene, seq))
									name = taxon+"_"+acc

									#if georeferencing required, check first
									if params.georef:
										found=False
										#print("checking for coordinates")
										for qual in gb_quals:
											if qual["GBQualifier_name"] == "lat_lon":
												if name and gene:
													name = name + "[" + qual["GBQualifier_value"] + "]"
												found=True
												continue
										if not found:
											continue

									#Add found to foundGenes
									if name and gene:
										if gene in foundGenes:
												foundGenes[gene][name] = seq
										else:
											foundGenes[gene] = dict()
											foundGenes[gene][name] = seq
									continue
						except Exception as e:
							if "GBSeq_primary-accession" in rec:
								l = list()
								if rec["GBSeq_organism"]:
									l.append(rec["GBSeq_organism"])
								else:
									l.append("NULL")
								if rec["GBSeq_definition"]:
									l.append(rec["GBSeq_definition"])
								else:
									l.append("NULL")
								bad_accessions[rec["GBSeq_primary-accession"]] = l
								print("Sequence",rec["GBSeq_primary-accession"],"lacks the proper XML fields! Skipping it and writing to file <bad_gis.log>:",e)

							else:
								print("Problem parsing record <NO READABLE ACCESSION>:",e)
				count += 1

				#Wait for pre-defined time before attempting next query
				time.sleep(params.time)


		if len(bad_accessions) > 0: #write log
			oname = params.out + "bad_gis.log"
			print("Writing bad accessions to",oname)
			list2file(bad_accessions, oname)
			#dictOlists2file(bad_accessions, "bad_gis.log")

		#if user only wants one sequene per species:
		keepOnes = dict()
		if params.keepOne:
			accessions = dict()
			#for each gene
			for gene in foundGenes:
				#find best sequence for each taxon
				for ind in foundGenes[gene]:
					#print(ind)
					spec = "_".join(ind.split("_")[:-1])
					acc = ind.split("_")[-1]
					#print(spec)
					#print(acc)
					if gene in keepOnes:
						if spec in keepOnes[gene]:
							#print(foundGenes[gene][ind])
							if len(keepOnes[gene][spec]) < len (foundGenes[gene][ind]):
								keepOnes[gene][spec] = foundGenes[gene][ind]
								accessions[gene][spec] = acc
						else:
							keepOnes[gene][spec] = foundGenes[gene][ind]
							accessions[gene][spec] = acc
					else:
						keepOnes[gene] = dict()
						accessions[gene] = dict()
						keepOnes[gene][spec] = foundGenes[gene][ind]
						accessions[gene][spec] = acc
		oname = params.out + "kept_accessions.log"
		print("Writing accessions of retained sequences to:",oname)
		if len(accessions) > 0: #write log
			d = pd.DataFrame(accessions)
			d.to_csv(oname, sep="\t", index=False)

		df = pd.DataFrame(foundGenes)
		#if only keeping one per species, replace DF with the selected bests
		if params.keepOne:
			df = pd.DataFrame(keepOnes)

		gnum = df.shape[1]
		tnum = df.shape[0]

		print("Found %s genes, %s taxa.\n" %(gnum, tnum))

		if params.taxCov:
			print("Removing taxa with less than %s data...\n" %params.taxCov)
			t = int(df.shape[1]*params.taxCov)
			if t < 1:
				t = 1
			df.dropna(thresh=t,how='all',axis=0,inplace=True)
		else:
			print("Dropping taxa with no data...\n")
			df.dropna(thresh=1,how='all',axis=0,inplace=True)

		if params.genCov:
			print("Removing genes with less than %s data...\n" %params.genCov)
			t = int(df.shape[0]*params.genCov)
			if t < 1:
				t = 1
			df.dropna(thresh=t,how='all',axis=1,inplace=True)

		if params.keepBest:
			print("Only retaining gene with most individuals...")
			gene_counts = dict()
			counts=df.count()
			print("Keeping gene:",counts.idxmax(), "\n")
			df.drop(df.columns.difference([counts.idxmax()]), 1, inplace=True)
			#print(df.count())


		print("Genes sampled:")
		if (df.shape[1] < 1 or df.shape[1] < 1):
			print("No data remains after filtering! Try lowering <-g> or <-t>!")
			sys.exit(0)
		print("Gene\t#taxa\t%Missing")
		for col in df:
			ct = df[col].count()
			t = len(df[col])
			print(col, "\t", ct , "\t", ((t-ct)/t)*100)

		#Write each gene alignment to a FASTA file
		print("\nPreparing to create gene alignments...")
		aln_files = list()
		for name, values in df.iteritems():
			outfile = params.out + name + ".fasta"
			aln_files.append(outfile)
			print("...Writing",outfile,"...")
			with open(outfile, 'w') as fh:
				try:
					for i,v in values.iteritems():
						if v not in ["nan", "NaN"] and pd.notnull(v):
							ol = ">" + i + "\n" + v + "\n"
							fh.write(ol)
				except IOError:
					print("Could not read file ",f)
					sys.exit(1)
				finally:
					fh.close()

		if params.noaln:
			print("\nSkipping alignment steps <-a,--noaln>")
			print("\nDone!\n")
		else:
			print("\nCalling MUSCLE on each alignment...")
			alignments = list()
			for aln in aln_files:
				name = aln.split(".fasta")[0]
				muscle_cline = MuscleCommandline(params.muscle, input=aln)
				stdout, stderr = muscle_cline()
				align = AlignIO.read(StringIO(stdout), "fasta")
				#trim leading and trailing gaps.
				# l = align.get_alignment_length()
				# print(name, "- Alignment length before trimming end gaps:",l)
				# #first pass through alignment to get characters to trim
				# lead = 0
				# trail = 0
				# for record in align:
				# 	rec_lead = countLeading(record.seq, "-")
				# 	rec_trail = countTrailing(record.seq, "-")
				# 	if rec_lead > lead:
				# 		lead = rec_lead
				# 	if rec_trail > trail:
				# 		trail = rec_trail
				#
				# align = align[:,lead:l-trail]
				# alignments.append(align)
				#
				# print("After removing leading and trailing gaps:",align.get_alignment_length())

				#Write NEXUS for alignment
				nex = name + "_muscle.nex"
				print("Writing alignment (nexus):",nex)
				align2nex(align, nex)

				fas = name + "_muscle.fasta"
				print("Writing alignment (fasta):",fas)
				align2fasta(align, fas)

			#Write NEXUS concatenated alignment

#function to write an AlignIO object as NEXUS
def align2nex(aln, nex):
	with open(nex, 'w') as fh:
		try:
			slen = aln.get_alignment_length()
			header = "#NEXUS\n\nBegin data;\nDimensions ntax=" + str(len(aln)) + " nchar=" + str(slen) + ";\n"
			header = header + "Format datatype=dna symbols=\"012\" missing=? gap=-;\nMatrix\n\n"
			fh.write(header)
			for rec in aln:
				sline = str(rec.id) + " " +str(rec.seq) + "\n"
				fh.write(sline)
			last = ";\nEnd;\n"
			fh.write(last)
		except IOError as e:
			print("Could not read file:",e)
			sys.exit(1)
		except Exception as e:
			print("Unexpected error:",e)
			sys.exit(1)
		finally:
			fh.close()

#function to write an AlignIO object as FASTA
def align2fasta(aln, nex):
	with open(nex, 'w') as fh:
		try:
			for rec in aln:
				line = ">" + str(rec.id) + "\n" +str(rec.seq) + "\n"
				fh.write(line)
		except IOError as e:
			print("Could not read file:",e)
			sys.exit(1)
		except Exception as e:
			print("Unexpected error:",e)
			sys.exit(1)
		finally:
			fh.close()


#Function to cound number of a leading character in string
def countLeading(s, c):
	return(abs(len(s.lstrip(c)) - len(s)))

#Function to cound number of a trailing character in string
def countTrailing(s, c):
	return(abs(len(s.rstrip(c)) - len(s)))

#function to write a list to file
def list2file(l, f):
	with open(f, 'w') as fh:
		try:
			for i in l:
				line = i + "\n"
				fh.write(line)
		except IOError:
			print("Could not read file ",f)
			sys.exit(1)
		finally:
			fh.close()

#function to write dict of lists to file
def dictOlists2file(d, outfile):
	with open(outfile, 'w') as fh:
		try:
			l1 = "GI\tTAXON\tDEFINITION\n"
			fh.write(l1)
			for key in d.keys():
				line = str(key) + "\t"
				for i in d[key]:
					# if i is None:
					# 	i="nan"
					line += i + "\t"
				line += "\n"
				fh.write(line)
		except IOError:
			print("Could not read file ",outfile)
			sys.exit(1)
		finally:
			fh.close()


#function to write dict of lists to file
def dictOdicst2file(d, outfile):
	with open(outfile, 'w') as fh:
		try:
			l1 = "GI\tTAXON\tDEFINITION\n"
			fh.write(l1)
			for key in d.keys():
				line = str(key) + "\t"
				for i in d[key]:
					# if i is None:
					# 	i="nan"
					line += i + "\t"
				line += "\n"
				fh.write(line)
		except IOError:
			print("Could not read file ",outfile)
			sys.exit(1)
		finally:
			fh.close()

#Fetch list of hits and write straight to file
def fetchToFile(hits, db, size, out, ft):
	outfile = out
	if ft == "fasta":
		print("Writing results to FASTA file...")
		outfile = out + ".fasta"
	else:
		print("Writing results to GENBANK file...")
		outfile = out + ".gb"

	with open(outfile, 'w') as OFH:
		try:
			#post NCBI query
			search_handle     = Entrez.epost(db, id=",".join(hits))
			search_results    = Entrez.read(search_handle)
			webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"]

			#fecth all results in batches
			for start in range(0,len(hits), size):
				handle = Entrez.efetch(db=db,retmode="text",rettype=type,retstart=start,retmax=size,webenv=webenv,query_key=query_key )
				OFH.write(handle.read())#get list of entries for given query
		except IOError:
			print("Could not read file ",outfile)
			sys.exit(1)
		finally:
			OFH.close()


#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'A:e:q:d:m:b:GFg:t:o:hm:aT:s', \
			["email=","query=","db=","retmax=","batch=","GBonly","Fonly",
			"genCov=", "taxCov=", "out=", "help", "muscle=", "noaln", "time=",'qlist=',
			"single", "georef", "bySpecies", "keepBest", "api="])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.email=None
		self.api=None
		self.query=None
		self.qlist=None


		self.db="nuccore"
		self.retmax=1000000000
		self.rettype='gb'
		self.batchSize=500
		self.out=""

		self.genCov=None
		self.taxCov=None


		self.GBonly = False
		self.Fonly = False

		self.muscle = "muscle"
		self.noaln = False

		self.keepOne=False
		self.georef=False
		self.perSpec=False
		self.keepBest=False

		self.time = 0.5

		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg in options:
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 'e' or opt == 'email':
				self.email = arg
			elif opt == 'h' or opt == 'help':
				pass
			elif opt == 'o' or opt == 'out':
				self.out = arg
			elif opt == 'q' or opt == 'query':
				self.query = arg
			elif opt == 'm' or opt == 'retmax':
				self.retmax = int(arg)
			elif opt == 'b' or opt == 'batch':
				self.batchSize = int(arg)
			elif opt == 'G' or opt.lower() == 'gbonly':
				self.GBonly = True
			elif opt == 'F' or opt.lower() == 'fonly':
				self.Fonly = True
			elif opt == 'g' or opt.lower() == 'gencov':
				self.genCov = float(arg)
			elif opt == 't' or opt.lower() == 'taxcov':
				self.taxCov = float(arg)
			elif opt == 'a' or opt == 'noaln':
				self.noaln = True
			elif opt == 'm' or opt == 'muscle':
				self.muscle = arg
			elif opt == 'T' or opt == 'time':
				self.time = float(arg)
			elif  opt == 'qlist':
				self.qlist = arg
			elif opt.lower() == "byspecies":
				self.perSpec=True
			elif opt.lower() == "georef":
				self.georef=True
			elif opt.lower() == "keepbest":
				self.keepBest=True
			elif opt == "s" or opt.lower() == "single":
				self.keepOne=True
			elif opt == "A" or "api":
				self.api=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.email:
			self.display_help("Email not provided <-e,--email>")
		if not self.query and not self.qlist:
			self.display_help("No query given! <-q,--query> or <--qlist>")

		#Any changes necessary
		if self.Fonly:
			self.retmode = 'fasta'

		if self.perSpec and self.keepOne:
			self.display_help("--single and --bySpecies cannot be used together.")

		if self.keepOne and self.georef:
			self.display_help("--single and --georef cannot be used together.")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfetcher.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: GenBank parser, fetches and aligns data for automated phylogenetic/phylogeography pipelines")
		print("""
	Mandatory arguments:
		-e, --email	: Email used for Entrez
		-A, --api	: API key (needed for >3 queries per second to avoid HTTP timeout)
		-q, --query	: Query encased in quotes (e.g. 'viperidae[Organism]')
		   -or-
		--qlist	: .txt file with a list of queries to perform

	Search arguments:
		-d,--db		: NCBI database to query (default: nuccore)
		-t,--time	: Time (seconds) to wait between batches (default: 0.5)
		-m,--retmax	: Maximum hits to return [default=1 bil.]
		-b,--batch	: Batch size for Entrex queries [default=500]
		-G,--GBonly	: ONLY fetch GenBank formatted hits and write to file
		-F,--Fonly	: ONLY fetch FASTA-formatted hits and write to file

	Data curation options:
		-s,--single	: Keep a single representative for each species
			By default, this is the longest (most complete) sequence
		-g,--genCov	: Proportion of taxa to keep a gene [default=OFF]
		-t,--taxCov	: Proportion of genes to keep a taxon [default=OFF]
			--Taxa must have a minimum of 1 sampled gene...
		--georef	: ONLY retain samples with coordinate data
			WARNING: Coordinate data won't be standardized. Proceed w/ caution
		--keepBest	: Only keep the gene with the best coverage per species

	Output options
		-o,--out	: Output prefix
		-m,--muscle	: Path to MUSCLE binary to align outputs (default=muscle)
		-a,--noaln	: Skip alignment step
		--bySpecies	: Output alignments for each species
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
