#!/soft/bin/python

import os,sys
import regex
from Bio import SeqIO
import scipy, scipy.cluster
import numpy as np
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from itertools import islice
import time
import subprocess
import glob
import operator

##This script should scan for all possible kmers and match it \
## with the control set.
## perl -one liner to get all possible k-mers of length k:
##perl -se '$,="\n"; print split / /, qx/echo @{["{$a}" x $k]}/;' -- -k=7 -a='A,T,G,C' >tmp

def main(argv):

	if len(sys.argv) < 2:
		sys.exit('Usage: script.py <cancer_type>')

	cancer = sys.argv[1]  #terminal input
	cancer = [cancer.lower()]

	#ev_type = ("a3","a5","mx","ri","se") #"al","af"
	ev_type = ["se"]
	#cancer = ("brca","coad","hnsc","kich","kirc","kirp","lihc","lusc","luad","prad","thca")
	regulation = ("positive", "negative")

	#define tools paths
	meme2images_path = "/soft/bio/sequence/meme-4.10.1_4/bin/meme2images"
	iupac2meme_path = "/soft/bio/sequence/meme-4.10.1_4/bin/iupac2meme"
	mosea_path="/projects_rg/babita/TCGA/soft/mosea/mosea.py"
	fimo_path = "/projects_rg/babita/TCGA/soft/meme/bin/fimo"

	#define dir paths
	tcga_dir = "/projects_rg/babita/TCGA/pipeline/run11"
	kmer_outdir = "/projects_rg/babita/TCGA/pipeline/run12/kmers_analysis"
	zscore_all_dir = "/projects_rg/babita/TCGA/pipeline/run12/kmers_analysis/zscores_all"
	
	if not os.path.exists(zscore_all_dir):
		os.makedirs(zscore_all_dir)

	#define parameters
	len_kmers = 7       #k-mers length to scan
	min_seq_cutoff = 90 #minimum percentile required for k-mers occurence in the sequences

	for can in cancer:
		for ev in ev_type:
			for reg in regulation:
				
				print("\n\nSTART Kmer-analysis for: {} | {} | {}".format(can,ev,reg))
				#create outdir per cancer type
				main_out_dir = os.path.join(kmer_outdir, can.lower(), ev.lower(), reg.lower())

				##scan the PFM to fasta sequqnce regulated & control
				scan_dir = os.path.join(main_out_dir, "scan") 
				enr_dir = os.path.join(main_out_dir, "enrichment") 

				fname = can.lower() + "_" + reg + ".fa"
				reg_fa_file =  os.path.join(tcga_dir, ev.upper(), "dpsi_files", fname)
				fname = fname + ".tab"
				reg_count_file = os.path.join(scan_dir, fname)
				

				fname = can.lower() + "_" + "control.fa"
				bg_fa_file = os.path.join(tcga_dir, ev.upper(), "dpsi_files", fname)
				fname = fname + ".tab"
				bg_count_file = os.path.join(scan_dir, fname)

				
   				fname = can.lower() + "_kmers_count.txt"
				kmers_count_outfile = os.path.join(main_out_dir,fname)


				fname = can.lower() + "_" + ev.lower() + "_" + reg.lower() + "_kmers_zscore.tab" 
				zscore_output_file = os.path.join(enr_dir, fname)

				if not os.path.exists(main_out_dir):
					os.makedirs(main_out_dir)

				if not os.path.exists(scan_dir):
					os.makedirs(scan_dir)

				if not os.path.exists(enr_dir):
					os.makedirs(enr_dir)


				'''
				call functions
				'''
				print("\n1.Function call to count Kmers on regulated set")
				count_kmers_on_sequence(reg_fa_file, kmers_count_outfile, len_kmers)

				print("\n2.Function call to create kmers pfms on regulated set")
				pfm_dir = analyse_kmers(main_out_dir, kmers_count_outfile, min_seq_cutoff,
					     meme2images_path, iupac2meme_path)

				print("\n3.Function call to scan motif and create pfms")
				motif_scan(reg_fa_file, pfm_dir, scan_dir, mosea_path, fimo_path)
				motif_scan(bg_fa_file, pfm_dir, scan_dir, mosea_path, fimo_path)

				print("\n4.Function call for enrichment analysis")
				enrichment_analysis(reg_fa_file, reg_count_file, 
           			bg_fa_file, bg_count_file, zscore_output_file, zscore_all_dir,
           			mosea_path)

				os.system("cp {} {}".format(output_file, zscore_all_dir))
				print("Done: output saved in: {}\n\n\n".format(zscore_output_file))


def _status_bar(count, count_all):
	a = (float(count)/count_all)*100
	a = int(round(a,2))
	b = int(round(float(a)/2, 0)) #determines the size of bar by printing '='

	sys.stdout.write('\r')
	sys.stdout.write("{}/{}".format(count,count_all) + "[%-50s] %d%%" % ('='*b, a) + "\t")
	sys.stdout.flush()


#k-mers sliding window
def slide_windows(seq, leng):
	seq = str(seq)
	seq = seq.rstrip()
	seq = seq.replace("--", "")
	it = iter(seq)

	result = tuple(islice(it, leng))

	if len(result) == leng:
		yield result

	for elem in it:
		result = result[1:] + (elem,)
		yield result
	


def count_kmers_on_sequence(reg_file, kmers_count_outfile, len_kmers):

	reg_fa_seq  = []
	reg_fa_id = []

	matched_motifs_sofar = []
	kmers_count = {}
	wf = open(kmers_count_outfile, 'w')

 
	#read regulated files matched for regions:
	for seq_record in SeqIO.parse(reg_file, "fasta"):
		reg_fa_id.append(seq_record.id)
		reg_fa_seq.append(seq_record.seq)

	
	zipped_reg_seq = zip(reg_fa_id, reg_fa_seq)

	#length of all regulated sequences
	fa_len = len(reg_fa_seq)


	print("\n\nScanning fasta sequences for kmers")
	cnt = 0
	for seq in reg_fa_seq:
		cnt += 1
		_status_bar(cnt, fa_len)

		kmer_list_per_sequence = []

		kmer_list_per_sequence = (["".join(x) for x in slide_windows(seq, len_kmers)])

		for kmer_found in kmer_list_per_sequence:

			if len(kmer_found) == len_kmers:

				if kmers_count.has_key(kmer_found):
					kmers_count[kmer_found] += 1
				else:
					kmers_count[kmer_found] = 1
	
	#sort dictonary created above on values (kmer counts)
	sorted_kmers = sorted(kmers_count.items(), key=operator.itemgetter(1))
	#print sorted_kmers

	for kmer, count in sorted_kmers:
		wf.write("{}\t{}\n".format(kmer, count))
	
	wf.close()
	print("\nkmers-count file created:\t{}\n\n".format(kmers_count_outfile))


def motif_scan(fa_file, pfm_path, out_dir, mosea_path, fimo_path):

	tmpdir = os.path.join(out_dir,"tmp2")

	os.system("/soft/devel/python-2.7/bin/python {} scan --scan --count --fasta {} \
		 	   --pfm_path {} --out_dir {} --fmo_path {} ".format(mosea_path,
		 	    fa_file, pfm_path, tmpdir, fimo_path))
	os.system("cp {}/* {}".format(tmpdir, out_dir))
	os.system("rm -rf {}".format(tmpdir))


def enrichment_analysis(reg_fa_file, reg_count_file,  bg_fa_file, bg_count_file,
						output_file,zscore_all_dir, mosea_path):

	os.system("/soft/devel/python-2.7/bin/python {} enrich --reg_fa_file {} \
		--reg_count_file {} --bg_fa_file {} --bg_count_file {} --out_file {}".format(
			mosea_path, reg_fa_file, reg_count_file, bg_fa_file, bg_count_file,
			output_file))


def analyse_kmers(main_dir, kmers_count_outfile, min_seq_cutoff,
				  meme2images_path, iupac2meme_path): 
	
	out_dir_logos = os.path.join(main_dir , "logos")
	out_dir_pfms = os.path.join(main_dir , "pfms")

	#if directory for logos and pfms & tmp doesnot exist, create
	if not os.path.exists(out_dir_logos):
		os.makedirs(out_dir_logos)

	if not os.path.exists(out_dir_pfms):
		os.makedirs(out_dir_pfms)


	#find minimum kmers cutoff
	#num_lines = sum(1 for line in open(kmers_count_outfile))
	#min_seq_cutoff = round((min_seq_cutoff * num_lines)/100)

	#find minimum quantile cutoff
	kmer_count=[]
	with open(kmers_count_outfile, 'r') as kcf:
		for line in kcf:
			s = line.split('\t')
			kmer_count.append(int(s[1]))

	min_cutoff_percentile = np.percentile(kmer_count, min_seq_cutoff)
	print("Input cutoff {} percentile.\nMinimum occurence selected {}\nMaximum occurence observed:{}".format(
		min_seq_cutoff,min_cutoff_percentile,max(kmer_count)))
	kcf.close()
	'''
	change date: 26/october/2015 
	Instead of clustering similar-kmers together, do enrichment analysis for each kmer 
	ocurring above the cut-off threshold.
	'''

	with open(kmers_count_outfile, 'r') as kcf:
		cnt = 1
		cnt2 = cnt
		for line in kcf:
			_status_bar(cnt2, len(kmer_count))
			s = line.split('\t')
			
			if int(s[1]) >= min_cutoff_percentile: #already considering min_frac_kmers option
				#seq.append(s[0])
				kmer = str(s[0]).rstrip()
				count = s[1].rstrip('\n')

				iupac = kmer

				pfm_file = iupac + "_" + str(count) + ".pfm" 
				pfm_file_path = os.path.join(out_dir_pfms, pfm_file)
				iupac_count = iupac + "_" + str(count)
			
				#4.create iupac id
				os.system("{} -numseqs {} {} >{} 2>/dev/null".format(iupac2meme_path, 
						count, iupac, pfm_file_path))
	
				#5.create Motif logo
				os.system("{} -motif {} -png {} {} 2>/dev/null".format(meme2images_path,
						iupac, pfm_file_path, out_dir_logos))

				cnt += 1
			cnt2 += 1

	print("\nPFMs and logo file created for : {} kmers".format(cnt))
	print("PFMs location: {}".format(out_dir_pfms))
	print("Logos location: {}\n".format(out_dir_logos))

	return out_dir_pfms
			
#To do :
#compare pfm files and seperate known and unknown motifs
#Do the enrichment test with controls
#repeat the analysis with pooled data


if __name__ == "__main__":
	main(sys.argv)