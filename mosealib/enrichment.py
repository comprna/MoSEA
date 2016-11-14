import sys
import os
import re
import fnmatch
from os.path import isfile, join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import itertools
import csv
from collections import defaultdict, deque
import math
import random
import numpy as np


'''
This script should read fasta file and motif count file and create a dictonaries, 
perform enrichment analysis and output zscores for each motif
'''
def main():
	#####################
	# required parameters
	reg_fa_file = "~/tmp/brca_se_positive_up.fa"
	reg_motif_count_file = "~/tmp/brca_se_positive_up.tab"
	bg_fa_file = "~/tmp/brca_se_control_up.fa"
	bg_motif_count_file = "~/tmp/brca_se_control_up.tab"
	output_file =  "~/tmp/output_motif_enr.txt"
	times_random = int(100)
	#length_correction_option = 0
	#length_percentage_correction = 15

	#call enrichment analysis
	calculate_enrichment(reg_fa_file, reg_motif_count_file,
						bg_fa_file, bg_motif_count_file,
						output_file, times_random )

	#print output
	print("Done enrichment, file saved in: {}".format(output_file))


def _check_file(filename):
	try:
		with open(filename): pass
	except IOError:
		print('Cannot open the file:'+ filename)
		sys.exit(1)


def _status_bar(motif_id, count, count_all):
	a = (float(count)/count_all)*100
	a = int(round(a,2))
	b = int(round(float(a)/2, 0)) #determines the size of bar by printing '='

	sys.stdout.write('\r')
	sys.stdout.write("{} {}/{}".format(motif_id,count,count_all) + "[%-50s] %d%%" % ('='*b, a))
	sys.stdout.flush()


def get_fasta_seq_dictonary(fa_file):
	#returns fasta files dictonary for length and gc content
	dict_fa_count = {}
	dict_fa_length = {}

	for seq_record in SeqIO.parse(fa_file, "fasta"):
		fa_id = seq_record.id
		faseq = seq_record.seq.upper()
		gc = GC(faseq)
		dict_fa_length[fa_id] = len(faseq) #length fasta seq
		#calculate gc contect distribution to nearest 10
		gc_content_decimal_distribution = math.floor(gc / 10) * 10 #10-bin window
		dict_fa_count[fa_id] = gc_content_decimal_distribution

	return dict_fa_length, dict_fa_count


def _get_motif_count_values_for_dict(motif_ids, fields, fa_count_dict,
									 fa_len_dict, count_file):
	#returns various features of sequence and motif counts for creating dictonary
	seq_id = fields.pop(0)
	length = int(fa_len_dict[seq_id])
	gc = float(fa_count_dict[seq_id])

	if(len(motif_ids) != len(fields)):
		print('Column headers did not match with number of fields in row {}'.format(seq_id))
		print('Please check file: {}'.format(fa_file))
		sys.exit()	
	
	motif_id_and_motif_count = zip(motif_ids, fields)

	return motif_id_and_motif_count, seq_id, length, gc


def get_motif_dictonary(motif_count_file, dict_fa_len, dict_fa_count):
	#create dictonary for motif
	dict_motif_count = []
	#read bg motif count file:
	motif_count =  open(motif_count_file, "r")
	r = csv.reader(motif_count, delimiter='\t')
	motif_ids = r.next()
	motif_ids.pop(0)

	for fields in r:
		#call function
		zipped_motif_id_and_motif_count, seq_id, length, gc = _get_motif_count_values_for_dict(
											   motif_ids, fields,												  
											   dict_fa_count,
											   dict_fa_len,
											   motif_count_file)

		#iterate through both motif id and motif count
		for motif_id, motif_count in zipped_motif_id_and_motif_count:
			dict_motif_count.append({"motif_id" : motif_id, "length" : int(length),
								 "gc" : int(gc), "motif_count": int(motif_count)}) 
	return dict_motif_count


def calculate_enrichment(reg_fa_file, reg_motif_count_file,
						bg_fa_file, bg_motif_count_file,
						output_file, times_random):

	#check files
	_check_file(reg_fa_file)
	_check_file(reg_motif_count_file)
	_check_file(bg_fa_file)
	_check_file(bg_motif_count_file)


	#output file:
	with open(output_file, 'w+') as wf:
		wf.write("rbp_id\tobserved\texpected.avg\texpected.sd\tz.score\n")

	'''
	Create dictonaries for control and regulated:
	1. fasta dictonary : seq length, GC count
	2. motif dictonary : motif_id, motif_count, seq_length, GC_count
	'''

	#get fasta files dictonary for control and regulated
	bg_dict_fa_length, bg_dict_gc_count = get_fasta_seq_dictonary(bg_fa_file)
	reg_dict_fa_length, reg_dict_gc_count = get_fasta_seq_dictonary(reg_fa_file)

	#create motif files dictonary for control and regulated
	background_motif_dict = get_motif_dictonary(bg_motif_count_file, 
							bg_dict_fa_length, bg_dict_gc_count)
	
	regulated_motif_dict = get_motif_dictonary(reg_motif_count_file, 
							reg_dict_fa_length, reg_dict_gc_count)

	#read reg motif count file:
	motif_count =  open(reg_motif_count_file , "r")
	r = csv.reader(motif_count, delimiter='\t')
	motif_ids = r.next()
	motif_ids.pop(0)

	cnt_all = len(motif_ids)
	cnt = 0

	print("Perform Motif enrichment for file: {}".format(reg_motif_count_file))
	
	for motif_id in motif_ids:
		
		cnt += 1
		_status_bar(motif_id, cnt, cnt_all)

		#outfile:
		wf = open(output_file, 'a+')
		
		result_random_all = []	
		result_observed_all = []
		motifs_binary_all = []

		
		for gc in range(10,100,10):
			#Matched regulated with background set (Expected)			
			test_motif_to_match = motif_id
			background_set = []
			regulated_set = []
			result_observed = []
			result_random = []		
			length_observed = 0

			background_set = [x for x in background_motif_dict if(
							 (x["gc"] == gc) and
							 (x["motif_id"]== motif_id)) ]

			regulated_set =  [x for x in regulated_motif_dict if(
							 (x["gc"] == gc) and
							 (x["motif_id"]== motif_id)) ]

			result_observed = [ x["motif_count"] for x in regulated_set]
			length_observed = len(result_observed)
			result_observed_all.extend(result_observed)
			

			if background_set == [] or regulated_set == []:
				next
			else:
				for i in range(0, length_observed):
					try:
						background_subset = np.random.choice(background_set, 
							times_random, replace=False)
					except ValueError:
						try:
							background_subset = np.random.choice(background_set, 
							times_random)
						except ValueError:
							background_subset=[]
							print("Not enough data found for gc_content:{}".format(
							gc))
							print("Skipping event")
							next

					motif_count_subset = []
					motif_count_subset = [ x["motif_count"] for x in background_subset]
					result_random.extend(motif_count_subset)

				#convert motif_count to only 'present' and 'absent'
				motifs_binary = []
				motifs_binary = [1 if i else 0 for i in result_random]
				motifs_binary_all.extend(motifs_binary)

			if motifs_binary_all != []:
				expected_count = float(sum(motifs_binary_all)/float(len(motifs_binary_all)))	
				result_random_all.append(expected_count)
			else:
				next
		
		#ignore Inf/NaN warnings
		np.seterr(divide='ignore', invalid='ignore')
		#find observed (regulated) motif present/ all events
		motifs_binary = []
		motifs_binary = [1 if i else 0 for i in result_observed_all]
		observed_count_total = sum(motifs_binary)/float(len(motifs_binary))	
		#expected (background or control) meand and sd
		#calculate average, sd, zscore of expected

		bg_expected_mean = np.mean(result_random_all)
		bg_expected_sd = np.std(result_random_all)

		#zscore = (observed - expected_mean) / expected_standard_deviation
		z_score = (observed_count_total - bg_expected_mean) / float(bg_expected_sd)

		#print("motif={},observed={},expectd.avg={},expectd.sd={}.zscore={}".format(
		#motif_id,observed_count_total,bg_expected_mean,bg_expected_sd,z_score))
		#
		#write to outfile
		wf.write("{}\t{}\t{}\t{}\t{}\n".format(motif_id,observed_count_total,
			bg_expected_mean,bg_expected_sd,z_score))
	
	wf.close()
		
	print("\nOutput enrichment file created: {}".format(output_file))


if __name__ == "__main__":
    main()
