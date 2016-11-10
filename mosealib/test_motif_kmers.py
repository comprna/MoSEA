#!/soft/bin/python


import regex
import operator
import difflib
import scipy, scipy.cluster
import numpy as np
from nltk import metrics, stem, tokenize
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, Motif
from Bio.Align.Applications import MuscleCommandline
import motility
from Bio.Alphabet import IUPAC, Gapped
from itertools import islice
import matplotlib.pyplot as plt
from matplotlib.pyplot import show
from hcluster import pdist, linkage, dendrogram
import time
from Bio.Align import AlignInfo
import subprocess
import glob

##This script should scan for all possible kmers and match it \
## with the control set.
## perl -one liner to get all possible k-mers of length k:
##perl -se '$,="\n"; print split / /, qx/echo @{["{$a}" x $k]}/;' -- -k=7 -a='A,T,G,C' >tmp

def main(argv):
	reg_fa_file = '/home/babita/kmer/kmer2/brca_reg_events.fa'
	bg_fa_file = '/home/babita/kmer/kmer2/brca_control_events.fa'
	out_dir = '/home/babita/kmer/kmer2/brca'
	kmer_file = '/home/babita/kmer/kmer2/seven-mers_motifs.txt'
	kmers_count_outfile = '/home/babita/kmer/kmer2/output_sorted_kmers'
	Rconsensus_script = "/projects_rg/babita/TCGA/soft/mosea/mosealib/iupac_concensus.R"

	muscle_path = "/soft/bio/sequence/muscle-3.8.31/muscle"
	meme2images_path = "/soft/bio/sequence/meme-4.10.1_4/bin/meme2images"
	iupac2meme_path = "/soft/bio/sequence/meme-4.10.1_4/bin/iupac2meme"


	len_kmers = 7       #k-mers length to scan
	min_seq_cutoff = 0.60 #minimum count required for k-mers presence in the sequences
	distance_cutoff = float(0.9) #minimum distance/mismatches allowed between two kmers to be clustered together 

	count_kmers_on_sequence(reg_fa_file, bg_fa_file,out_dir, kmer_file, 
							kmers_count_outfile, len_kmers)

	analyse_kmers(out_dir, kmers_count_outfile, 
				 	min_seq_cutoff, distance_cutoff, muscle_path,
				 	meme2images_path, iupac2meme_path, 
				 	Rconsensus_script)


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
	


def count_kmers_on_sequence(reg_file, bg_file, out_dir, kmer_file, 
							kmers_count_outfile, len_kmers):

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
	print("\nkmers-count file created:\t{}".format(kmers_count_outfile))


#Input: z= linkage matrix, treshold = the treshold to split, n=distance matrix size
def split_into_clusters(link_mat,thresh,n):  #http://www.blopig.com/blog/2014/05/distance-matrix-clustering/
	c_ts=n
	clusters={}
	for row in link_mat:
		#print row
		if row[2] < thresh:
			n_1=int(row[0])
			n_2=int(row[1])

			if n_1 >= n:
				link_1=clusters[n_1]
				del(clusters[n_1])
			else:
				link_1= [n_1]

			if n_2 >= n:
				link_2=clusters[n_2]
				del(clusters[n_2])
			else:
				link_2= [n_2]

			link_1.extend(link_2)
			clusters[c_ts] = link_1
			c_ts+=1
		else:
			return clusters

def enrichment_analysis(reg_fa_file, reg_count_file, 
           	bg_fa_file, bg_count_file, output_file):
	mosea_path="/projects_rg/babita/TCGA/soft/mosea/mosea.py"

	os.system("python {} enrich --reg_fa_file {} \
			--reg_count_file {} \
           --bg_fa_file {} \
           --bg_count_file bg_sample_output.space \
           --out_file {} \
           --match_len 1 \
           --len_ext 20".format(mosea_path, reg_fa_file, reg_count_file, 
           	bg_fa_file, bg_count_file, output_file))

def analyse_kmers(out_dir, kmers_count_outfile, min_seq_cutoff,
					distance_cutoff, muscle_path, meme2images_path, 
					iupac2meme_path, Rconsensus_script_path):

	if out_dir == '':
		out_dir = './'
	
	out_dir_logos = os.path.join(out_dir + "/logos")
	out_dir_pfms = os.path.join(out_dir + "/pfms")

	#if directory for logos and pfms & tmp doesnot exist, create
	if not os.path.exists(out_dir_logos):
		os.makedirs(out_dir_logos)

	if not os.path.exists(out_dir_pfms):
		os.makedirs(out_dir_pfms)

	if not os.path.exists("tmp_dir"):
		os.makedirs("tmp_dir")

	os.system('rm tmp_dir/*')
	os.system('rm {}/*'.format(out_dir_logos))	
	os.system('rm {}/*'.format(out_dir_pfms))
	
	seq=[]
	count=[]
	dict_kmer_counts = []
	
	#tutorial : http://codereview.stackexchange.com/questions/37026/string-matching-and-clustering
	#tutorial: longest substring: http://www.bogotobogo.com/python/python_solutions.php#longest_substring
	
	#find minimum kmers cutoff
	num_lines = sum(1 for line in open(kmers_count_outfile))
	min_seq_cutoff = round((min_seq_cutoff * num_lines)/100)

	'''
	change date: 26/october/2015 
	Instead of clustering similar-kmers together, do enrichment analysis for each kmer 
	ocurring above the cut-off threshold.
	'''
	#find clusters of matched kmers
	
	with open(kmers_count_outfile, 'r') as kcf:
		for line in kcf:
			s = line.split('\t')
			
			if int(s[1]) >= min_seq_cutoff: #already considering min_frac_kmers option
				#seq.append(s[0])
				kmer = s[0]
				count = s[1]
				#create dict with seq and count
				#dict_kmer_counts.append({"kmer" : kmer, "count" : int(count)})
				'''
				code instered on 26-october-2015
				for each kmer abover cut-off do the enrichment analysis
				'''
				#test files:
				reg_fa_file = "/home/babita/kmer/kmer2/brca_se_negative.fa"
				reg_count_file = "/home/babita/kmer/kmer2/brca_se_negative.space"
				bg_fa_file = "/home/babita/kmer/kmer2/brca_se_control.fa"
				bg_count_file = "/home/babita/kmer/kmer2/brca_se_control.tab"
				output_file = "/home/babita/kmer/kmer2/brca_se_negative_enrichment.tab"

				enrichment_analysis(reg_fa_file, reg_count_file, 
           			bg_fa_file, bg_count_file, output_file):

	'''
	#old code:		
	mat_len = len(seq)
	matrix = [[0 for x in range(mat_len)] for x in range(mat_len)]

	i = 0
	j = 0

	start = time.time()
	#print start

	#print("\n\nClustering kmers (minimum required: {})".format(min_seq_cutoff))
	for i in range(mat_len):
		ii= i+1
		_status_bar(ii, mat_len)

		for j in range(mat_len):
			s1 = seq[i]
			s2 = seq[j]
			seq_diff = difflib.SequenceMatcher(None,s1,s2)
			ratio = seq_diff.ratio()
			matrix[i][j] = ratio
	
	
	#cluster based on distance matrix
	
	z = linkage(matrix, method='ward')
	#z = linkage(matrix, method='weighted') #try different methods
	#z = scipy.cluster.hierarchy.fclusterdata(matrix,1.0)

	#distance_cutoff = minimum distance/mismatches allowed between two
	#kmers to be clustered together
	
	cutoff = distance_cutoff  
	rows = mat_len

	clustering = split_into_clusters(z, cutoff, rows)

	if clustering==None:
		print("No clusters! distance cutoff is too high")
		print("All kmers are in the same cluster!")
		quit()

	#print clustering
	#quit()
	motif_cluster = []

	cluster_num = 0
	#print("\n\nAnalyzing kmer clusters")
	
	for cluster in clustering:
		cluster_num+=1
		mot = []		
		_status_bar(cluster_num, len(clustering))

		for i in clustering[cluster]:
			mot.append(seq[i])		
		
		motif_cluster.append(mot)
	
	print motif_cluster

	For Each cluster (c) thus created in motif_cluster, do following:
		1-create a fasta file with all k-mers in cluster (c)
		2-perform alignment of this fasta file using muscle
		3-create IUPAC code for DNA using this alignment file 
		4-create PFMs for Meme using IUPAC code generated
		5-create motif logo using PFM file generated
	

	cnt = 0

	for sequences in motif_cluster:
		cnt += 1
		kcount = 0

		for s in sequences:
			kmer_count = [x for x in dict_kmer_counts if(
						x["kmer"] == s)]
			k_count = [x["count"] for x in kmer_count]
			kcount += int(k_count[0])
	'''


		#1.create a fasta file with all k-mers in cluster (c)
		
		records = (SeqRecord(Seq(s)) for s in sequences)
		SeqIO.write(records, "tmp_dir/out_kmers_cluster.fasta", "fasta")


		align_file = "aligned_kmers_cluster_" + str(kcount) + ".fa"
		align_file_path = os.path.join("tmp_dir", align_file)

		#2. perform alignment of this fasta file using muscle
		os.system("{} -in tmp_dir/out_kmers_cluster.fasta -out {} -refine 2>/dev/null".format(
			muscle_path, align_file_path))

		aligned_kmers = []
		for seq_record in SeqIO.parse(align_file_path, "fasta"):
			s = str(seq_record.seq)
			#s = s.replace("-", "N")
			aligned_kmers.append(s)
			
		#3.
		try:
			iupac = motility.make_iupac_motif(aligned_kmers)
			print("{}".format(aligned_kmers))

		except Exception, err:
			#for gapped alignment
			print err
			print("{}".format(aligned_kmers))
			output = subprocess.Popen(['Rscript', Rconsensus_script_path, 
					 align_file_path], shell=False, stdout=subprocess.PIPE,
					 stderr=subprocess.PIPE)
			iupac = output.stdout.read()

		pfm_file = iupac + "_" + str(kcount) + ".pfm" 
		pfm_file_path = os.path.join(out_dir_pfms, pfm_file)
		iupac_count = iupac + "_" + str(kcount)
		
		#4.
		os.system("{} -numseqs {} {} >{} 2>/dev/null".format(iupac2meme_path, 
					kcount, iupac, pfm_file_path))
		#5.
		os.system("{} -motif {} -png {} {} 2>/dev/null".format(meme2images_path,
					iupac, pfm_file_path, out_dir_logos))

		print("\nDone:{}\tPFMs and logo file created for : {}".format(cnt,
					iupac_count)) 

		#Delete temperory fasta file & alignment file (optional)
		os.system("rm tmp_dir/out_kmers_cluster.fasta")

	print("output saved in dir: {}".format(out_dir))
#To do :
#compare pfm files and seperate known and unknown motifs
#Do the enrichment test with controls
#repeat the analysis with pooled data


if __name__ == "__main__":
	main(sys.argv[1:])
