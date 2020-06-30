import sys
import os
import re
import fnmatch
import numpy as np
from os.path import isfile, join
from Bio import SeqIO
from Bio.SeqUtils import GC
from time import sleep
import pandas as pd

try:
	from StringIO import StringIO ## for Python 2
except ImportError:
	from io import StringIO ## for Python 3

try:
	import commands	## for Python 2
except ImportError:
	import subprocess ## for Python 3


def main():
	pass
	#print("In main")

def _status_bar(count, count_all):
	a = (float(count)/count_all)*100
	a = int(round(a,2))
	b = int(round(float(a)/2, 0)) #determines the size of bar by printing '='

	sys.stdout.write('\r')
	sys.stdout.write("{}/{}".format(count,count_all) + "[%-50s] %d%%" % ('='*b, a) + "\t")
	sys.stdout.flush()


def _check_file(filename):
	try:
		with open(filename): pass
	except IOError:
		print('Cannot open the file:'+ filename)
		sys.exit(1)

def _which(file):
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
                return path + "/" + file

    return None


def _create_motif_count_list(motif_count_file, df_seq):	

	_check_file(motif_count_file)

	#create dictonary for motif
	list_fmo_fa_count = []

	df = pd.read_table(motif_count_file, sep = "\t", header=0)
	df = df[['sequence name']]
	df.columns = ['event']

	df['count'] = df.groupby('event')['event'].transform(pd.Series.value_counts)
	df = df.drop_duplicates()

	df = pd.merge(df_seq, df,  on="event", how = 'outer')
	df = df.fillna(0)

	dflist_count = df['count'].tolist()
	return dflist_count

############	

def get_fasta_sequence(bedtools_path, bed_file, genome, out_file):
	'''
	to extract fasta sequences from given genome coordinates
	and genome.fasta file
	'''

	_check_file(bed_file)

	#out_file: getbase name, remove '.bed' write '.fa'
	bed_name = os.path.basename(bed_file)
	bed_name = re.sub('\.bed$', '', bed_name) #bed_name.replace('.bed','')
    
	os.system("{} getfasta -fi {} -bed {} -s -name -fo {}".format(
			 bedtools_path, genome, bed_file, out_file))



def scan_motif_pfm(fa_file, out_dir, pfm_dir, fimo_path, count_flag):
	'''
	to scan motif on fasta seuqnces using FIMO binary
	'''

	_check_file(fa_file)

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	if _which('fimo') is None and os.path.exists(fimo_path) is False:
		print('Script need Fimo from the MEME suite in your path, available at: http://meme.nbcr.net/meme/')
		print("Cannot run without fimo. Exiting..!!")
		sys.exit(1)

	status, fimo_version = commands.getstatusoutput("{} --version".format(fimo_path))

	if os.path.isdir(pfm_dir):
		pfmfiles = [f for f in os.listdir(pfm_dir) if isfile(os.path.join(pfm_dir,f)) ]
		cnt_all = len(pfmfiles)

	elif os.path.exists(pfm_dir):
		pfmfiles = [pfm_dir]
		cnt_all = 1
	else:
		print("PFM file/dir : {} NOT FOUND".format({pfm_dir}))
		quit()

	fa_name = os.path.basename(fa_file) #basename

	cnt = 0

	pval = '1e-3'  #Default p-value, optional to change

	print("scanning Motifs on file: {}".format(fa_file))

	for i in range(cnt_all):

		cnt += 1
		_status_bar(cnt, cnt_all)

		pfm = pfmfiles[i]

		motif_name = os.path.basename(pfm)
		#motif_name = os.path.splitext(pfm)[0] #basename

		motif_name = motif_name.split(".")[0]

		if os.path.isdir(pfm_dir):
			pfm = os.path.join(pfm_dir,pfm)

		fmo_file_name = fa_name + "." + motif_name + ".fmo"

		#if fimo_version == '4.8.1':
		try:
		##elif fimo_version in ['4.10.1','4.10.1_1','4.10.1_2','4.10.1_3','4.10.1_4',"4.9.0"]:
			os.system("{} --norc --text --thresh {} {} {} > {} 2>/dev/null".format(fimo_path,
				pval, pfm, fa_file, os.path.join(out_dir,fmo_file_name)))

		# except Exception, err:
		# #else:
		# 	print "\nError!! Not a valid version of FIMO (MEME-suite). Please install meme-4.9.0 or above\n\n"
		# 	print("Exception Error {}".format(err))
		# 	sys.exit()

		except Exception as error:
			print("\nError!! Not a valid version of FIMO (MEME-suite). Please install meme-4.9.0 or above\n\n")
			print("Exception Error {}".format(err))
			sys.exit()

	print("\nScanned {} motif(s). Output saved in dir: {}".format(cnt_all, out_dir))


	#if '--count' option is selected
	if count_flag:
		count_motif(fa_file, pfm_dir, out_dir)


def scan_motif_kmer(fa_file, out_dir, kmer_dir, count_flag):
	'''
	to scan motif on fasta seuqnces using kmers file(s)
	'''
	_check_file(fa_file)


	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	if os.path.isdir(kmer_dir):
		kmerfiles = [f for f in os.listdir(kmer_dir) if isfile(os.path.join(kmer_dir,f)) ]
		cnt_all = len(kmerfiles)

	elif os.path.exists(kmer_dir):
		kmerfiles = [kmer_dir]
		cnt_all = 1

	else:
		print("Kmer file or directory : {} NOT FOUND".format({kmer_dir}))
		quit()

	fa_name = os.path.basename(fa_file) #basename

	cnt = 0

	print("scanning Motifs on file: {}".format(fa_file))

	for i in range(cnt_all):

		cnt += 1
		_status_bar(cnt, cnt_all)

		kmerfile = kmerfiles[i]

		motif_name = os.path.basename(kmerfile)
		motif_name = motif_name.split(".")[0]

		if os.path.isdir(kmer_dir):
			kmerfile = os.path.join(kmer_dir, kmerfile)

		fmo_file_name = fa_name + "." + motif_name + ".fmo"
		fmo_file_name = os.path.join(out_dir, fmo_file_name)

		#open output file
		with open(fmo_file_name, 'w'): pass


		# Check if file is fasta file ot not
		flag = 0
		for line in open(kmerfile):

			if line.startswith(">"):
	        
				flag = 1
				break

		#==================================================================#
		kmer_seq = []

		if flag:
			#read fasta files of kmers
			kmer_fasta_seq = SeqIO.parse(open(kmerfile), "fasta")

			for seq_record in kmer_fasta_seq :
				kmer_seq.append(str(seq_record.seq))
		
		else :

			for line in open(kmerfile):
				li = line.strip()
				kmer_seq.append(li)

		#convert all lower case string to upper case
		kmer_seq = [f.upper() for f in kmer_seq]

		#replace U with T's if the kmers are RNAseq
		kmer_seq  = [f.replace('U', 'T') for f in kmer_seq]


		#===================================================================#
		#open outfile in append mode
		of = open(fmo_file_name,'a')
		of.write("#pattern name\tsequence name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched sequence\n")

		#For main Fasta sequqnces to scan
		fa_seq = []
		fasta_seq = SeqIO.parse(open(fa_file), "fasta")


		for seq_record in fasta_seq:
			
			fa_seq = str(seq_record.seq)
			fa_header = seq_record.id

			#convert all lower case string to upper case
			fa_seq = fa_seq.upper()
			fa_seq = fa_seq.replace('U', 'T') 
	

			#Scan kmer_seq on fa_seq and create output file in fmo format:
			for kmer in kmer_seq:
				Q = re.compile(kmer)
				found = [item.start(0)+1 for item in Q.finditer(fa_seq)]
				len_kmer = len(kmer)

				for f in found:				
					of.write("{}\t{}\t{}\t{}\t+\t0\t0\t\t{}\n".format(motif_name, fa_header, f, (f+len_kmer)-1, kmer))
				
		of.close()

	print("\n")
	#if '--count' option is selected, create count table
	if count_flag:
		count_motif(fa_file, kmer_dir, out_dir)


def count_motif(fa_file, pfm_dir, out_dir):
	
	'''
	Creates tab-delimited count file for all PFMs or Kmers of the motif provided in the folder
	'''
	print(fa_file, pfm_dir, out_dir)
	_check_file(fa_file)

	##create list of fasta sequqnces
	fa_seq_ids  = []
	for seq_record in SeqIO.parse(fa_file, "fasta"):
		fa_seq_ids.append(seq_record.id)

	
	df = pd.DataFrame({"event" : fa_seq_ids})
	rowlen = len(fa_seq_ids) 

	pfmfiles = []

	if os.path.isdir(pfm_dir):
		#Read all motif fmo files 
		pfmfiles = [f for f in os.listdir(pfm_dir) if (isfile(join(pfm_dir, f)) and os.path.getsize(join(pfm_dir,f)) ) ]
	
	elif os.path.exists(pfm_dir) and os.path.getsize(pfm_dir): #if it is just one file instead of a directory
		pfmfiles = [pfm_dir]	


	fa_name = os.path.basename(fa_file) #basename
	out_file = fa_file + ".tab"

	fmofiles = []


	for pfm in pfmfiles:

		if os.path.isdir(pfm_dir):
			pfm_name =  os.path.basename(pfm) ###edit the path,
			pfm_name = pfm_name.split(".")[-2]  #remove trailing ".pfm"
		else:
			pfm_name = pfm.split(".")[-2] #remove trailing ".pfm"
			

		scanned_fmo_file_name = fa_name + "." + pfm_name + ".fmo"

		#print(os.path.join(out_dir, scanned_fmo_file_name))

		if os.path.getsize(os.path.join(out_dir, scanned_fmo_file_name)):
			fmofiles.append(scanned_fmo_file_name)
	

	if not len(fmofiles):
		print("\n\nNone of the provided motifs found on the Fasta sequences. Exiting!\n\n")
		quit()
	

	fmofiles = sorted(fmofiles)
	collen = len(fmofiles)
	collen += 1 #adding 1 to no. of cols, for 'event'

	##read all fmo files one by one and create dict of all fmo motif count
	motif_count_list = []
	cnt = 0
	collen_org = collen - 1 #because 1 was added earlier, to original count

	print("Counting Motifs on file: {}".format(fa_file))

	for fmofile in fmofiles:
		
		cnt += 1
		_status_bar(cnt, collen_org)

		motif_id = fmofile.split(".")[-2]	
		fmofile = os.path.join(out_dir,fmofile)

		motif_count_list = _create_motif_count_list(fmofile, df)
	
			
		if len(motif_count_list) == len(fa_seq_ids):
			df[motif_id] = pd.Series(motif_count_list, index=df.index)
			
		else:
			print("\nLength mismatch error in creating count table:")
			print("Length of motif counts set = {}".format(len(motif_count_list)))
			print("Length of fasta sequence ids set = {}".format(len(fa_seq_ids)))
			print("Please check for motif_id : {}".format(motif_id))
			print("Exiting!!\n")
			sys.exit()

		
	if (df.empty or (df.shape[1] != collen) or (df.shape[0] != rowlen)):
	 	print('\nError in creating table. Count file not created!') 

	else:
		#if all is well, write the df into table:
		df.to_csv(out_file, sep='\t', index=False)
	
	print("\nFile created: {}".format(out_file))


if __name__ == "__main__":
    main()
