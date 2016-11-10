'''
#
# merge *.fmo output and MAF files from TCGA
# MAF files description: https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
# fmo files: output file gernerated from fimo (meme package) motif scan, converted to bed format by running script: fmo_to_bed_for_suppa_events
#fmo.bed format: Should have alteast six compulsory columns:
Compulsory Columns:
	chr7	77523366	77523372	PHTF2|57157;SE:chr7:77522981-77523299:77523310-77530028:+;down	200	+
	1.chromosome, 2.motif_start_position, 3.Motif_end_position, 4.seq_id, 5.length_of_seq, 6.strand
Optional columns:
	A1CF_00001	TTAATTG	77523310	77523510
	7.Motif_id, 8.Motif_kmer, 9.Sequence_start_postion, 10.Sequence_end_postion
'''


##read MAf files, seperate desired columns.
##merge fmobed and maf files.


import sys,os 
from os.path import isfile, join
import re


def main():

	input_len = len(sys.argv)

	if input_len < 3:
		print "Usage: script.py <input_file*.fmo.bed> <input_maf_file.maf> <out_dir(optional)>"
		print "Help:\n\tCreates merged fmobed & maf file from motif coordinates from *.fmo.bed and *.maf input file"
		print "\tOutput: fmo.bed.maf file (input_file_name.maf\n\n)"
		sys.exit()

	
	input_fmobed_dir = sys.argv[1]
	input_maf_file = sys.argv[2]
	output_file = sys.argv[3]


	#if input_len == 3:	
	#	output_file = sys.argv[3]
	#	if not os.path.exists(output_file):
	#		os.makedirs(output_file)	
	#else:
	#	output_file = 'merged_fmobed_maf_file.tab'


	'''
	Read input file, check no. of columns (should be in fimo format)
	'''
	#basename = os.path.basename(input_fmobed_dir) #basename of *.fmo.bed
	#output_file = basename + ".maf"

	basename_maf = os.path.basename(input_maf_file)
	cancer_type = basename_maf.split("_")[0]

	#if out_dir != '':
	#	output_file = os.path.join(out_dir, output_file)

	fo = open(output_file, "w")
	#add header
	headers = []
	headers = ("cancer","motif_id","event","kmer","seq_len","chromosome",
					"sequence_start","sequence_end","strand",
					"motif_start","motif_end","mutation_start",
					"mutation_end","variant_class","variant_type","ref_allele","tumor_allele_1",
					"tumor_allele_2","pateint","sample_barcode")

	fo.write('\t'.join(headers) + '\n')
	#print('\t'.join(headers))

	#if os.path.exists(output_file) is False:

		#print("Directory: {} doesn't exist".format(out_dir))
	#	print("Path for outfile: {} not found".format(output_file))
	#	sys.exit()
	

	maffields_dict = []

	with open(input_maf_file) as fi:
		for line in fi:
			line.strip()
			maffields = line.split('\t')
			##assign values to the field

			maffields_dict.append({"maf_gene" : maffields[0],
									"maf_entrez" : maffields[1],
									"maf_chr" : ("chr"+maffields[4]),
									"maf_chr_start" : int(maffields[5]),
									"maf_chr_end" : int(maffields[6]),
									"maf_chr_strand" : maffields[7],
									"maf_variant_class" : maffields[8],
									"maf_variant_type" : maffields[9],
									"maf_ref_allele" : maffields[10],
									"maf_tumor_allele1" : maffields[11],
									"maf_tumor_allele2" : maffields[12],
									"maf_sample_id" : maffields[15]
									})

	fi.close()


	fmobed_dict = []

	input_fmobed_files = [f for f in os.listdir(input_fmobed_dir) if (isfile(os.path.join(input_fmobed_dir,f)) and f[-8:] == ".fmo.bed") ]
	if len(input_fmobed_files) == 0:
		print("No fmo.bed files found. All fmobed files should have extension '*.fmo.bed' ")
		sys.exit()

	for input_fmobed_file in input_fmobed_files:
		input_fmobed_file = os.path.join(input_fmobed_dir, input_fmobed_file)
		with open(input_fmobed_file) as fi:	
			for fields in fi:
				field = fields.strip('\n')	
				match_fmobed_coordinates(field, maffields_dict, cancer_type, fo)

		fi.close()
	fo.close()

######
#match the motif coordinates from fmo bedfiles and write it in the merge file

def match_fmobed_coordinates(line, maf_fields_dict, can, fh):

	fmofields = line.split('\t')
	matched_set = []

	if len(fmofields) != 10:
		print("Not correct number of fields in fmobed file:")
		print('\t'.join(fmofields))
		sys.exit()
	
 	fmo_chr = fmofields[0]
	fmo_motif_start = int(fmofields[1])
	fmo_motif_end = int(fmofields[2])
	fmo_eventid = fmofields[3]
	fmo_len = fmofields[4]
	fmo_strand = fmofields[5]
	fmo_motif = fmofields[6]
	fmo_kmer = fmofields[7]
	fmo_org_start = fmofields[8]
	fmo_org_end = fmofields[9]


	matched_set = [x for x in maf_fields_dict if(
							 (x["maf_chr"] == fmo_chr) and
							 (x["maf_chr_strand"] == fmo_strand) and
							 (x["maf_chr_start"] >= fmo_motif_start) and
							 (x["maf_chr_end"] <= fmo_motif_end)) ]

	if len(matched_set):
		outlist_fmo = []
		outlist_fmo.extend((can,fmo_motif,fmo_eventid,fmo_kmer,fmo_len,fmo_chr,
			fmo_org_start,fmo_org_end,fmo_strand,str(fmo_motif_start),str(fmo_motif_end)))


		for x in matched_set:
			
			outlist_maf = []
			a = x["maf_chr_start"]
			b = x["maf_chr_end"]
			outlist_maf.append(str(a))
			outlist_maf.append(str(b))
			outlist_maf.append(x["maf_variant_class"])
			outlist_maf.append(x["maf_variant_type"])
			outlist_maf.append(x["maf_ref_allele"])
			outlist_maf.append(x["maf_tumor_allele1"])
			outlist_maf.append(x["maf_tumor_allele2"])
			outlist_maf.append(x["maf_sample_id"])
			a = x["maf_sample_id"]
			a = a.split('-')
			sample_id = a[2] + "T"
			outlist_maf.append(sample_id)

			fh.write('\t'.join(outlist_fmo) + '\t' )
			fh.write('\t'.join(outlist_maf) + '\n' )

if __name__ == "__main__":
    main()
