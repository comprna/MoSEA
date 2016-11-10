import sys, os
'''
FMO Extractor:
	Converts *.fmo files to bed files with motif coordinates.
	input format (fimo):
	#pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence
	A1CF_00001	AACSL|729522;A3:chr5:178195668-178199430:178195665-178199430:-;V2	163	169	+	8.30709	0.000862		CTAATTA
	output format (bed):
	#chr 	start_pos 	end_pos		event_id;motif_id;motif_kmer	strand
'''


#1. get input file
def main():

	input_len = len(sys.argv)

	if input_len < 2:
		print "Usage: script.py <input_fmo_file*> <out_dir(optional)>"
		print "Help:\n\tCreates bed file for motif coordinates from *.fmo input file"
		print "\tOutput: bed file (input_file_name.bed)"
		sys.exit()

	elif input_len > 2:	
			input_file = sys.argv[1]
			out_dir = sys.argv[2]	
	else:
		input_file = sys.argv[1]
		out_dir = ''

	'''
	Read input file, check no. of columns (should be in fimo format)
	'''

	basename = os.path.basename(input_file) #basename of *.fmo
	output_file = basename + ".bed"

	if out_dir != '':
		output_file = os.path.join(out_dir, output_file)

	fo = open(output_file, "w")
	fo.write("")

	if os.path.exists(output_file) is False:
		print("Directory: {} doesn't exist".format(out_dir))
		print("Path for outfile: {} not found".format(output_file))
		sys.exit()

	
	with open(input_file) as fi:
		for line in fi:
			line = line.strip()
			if line.startswith("#") or line.startswith("Motif"): #header ignore
				continue
					
			'''
			count no. of columns, check for event_id in second column
			'''
			try:
				line1 = line.split('\t')
				len_arr = len(line1)
				last_elem = len_arr -1
				motif_id = line1[0]
				event = line1[1]
				start = line1[2]
				end =  line1[3]
				kmer = line1[last_elem]

				'''
				get event_type (A3) from event_id by spliting ';'..':'
				example:AACSL|729522;A3:chr5:178195668-178199430:178195665-178199430:-;V2
				'''
				#print event
				ev = event.split(';')[1]
				event_type = ev. split(':')[0]

				'''
				get variable region (e.g V2) from last element of the event_id.
				'''
				var_region = event.rsplit(';', 1)[1]
				'''
				pass the values to the functions based on event_type
				'''

				call_event_func = {
					"A3" : fun_A3_bedfile,
					"A5" : fun_A5_bedfile,
					"SE" : fun_SE_bedfile,
					"MX" : fun_MX_bedfile,
					"RI" : fun_RI_bedfile,
					"AF" : fun_AF_bedfile,
					"AL" : fun_AL_bedfile
					}

				call_event_func[event_type.upper()](event, event_type, var_region,
													motif_id, int(start), int(end), kmer,
													output_file)				

			except Exception, err:
				print("Error occurred while parsing input file:\t{}".format(err))
				sys.exit()
		fi.close()
		fo.close()



def fun_A3_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	##ref: see Fig.2 suppa documentaion of A3 e1 & s1..descriptions
	#https://bitbucket.org/regulatorygenomicsupf/suppa
	#Example event id: (e1-s2:e1-s3)
	#ABCG2|9429;A3:chr4:89016761-89018605:89016751-89018605:-
	-ve strand: (s2-e1:s3-e1)
	'''
	fo = open(out_file, "a")
	len_ext = 200

	ev_all = event.split(';')	
	ev_coord = ev_all[1]
	ev_all = ev_coord.split(':')
	ev_type = ev_all[0]
	ev_chr = ev_all[1]
	ev_strand = ev_all[4]

	'''
	double check
	ev_type should be equal to event_type passed from the main function
	'''
	if event_type != ev_type:
		print("Error in event_id: {}, event type {} and {} did not match".format(
				event, event_type, ev_type))
		sys.exit()

	if(ev_strand=="-"):
		#eg.ABCG2|9429;A3:chr4:89016761-89018605:89016751-89018605:-   
		ev_s2_e1 = ev_all[2]
		ev_s2_e1 = ev_s2_e1.split('-')
		s2 = int(ev_s2_e1[0])	
		e1 = int(ev_s2_e1[1])	

		ev_s3_e1 = ev_all[3]
		ev_s3_e1 = ev_s3_e1.split('-')
		s3 = int(ev_s3_e1[0])

		if var_region.upper() == 'V1':
			e = len_ext - start 
			s = len_ext - end
			motif_coord1 = s3 + (s+1) #one-base correction previously done
			motif_coord2 = s3 + (e+1)	#,,
			org_start =  s3  
			org_end = s3 + len_ext

		elif var_region.upper() == 'V2':
			motif_coord1 = (s3 - end) + 1
			motif_coord2 = (s3 - start) + 1
			org_start =  s3 - len_ext
			org_end = s3

		else:
			print("var_region {} not defined. Error!".format(var_region))
			quit()

	elif(ev_strand=="+"):  
		#eg. A2BP1|54715;A3:chr16:7568391-7629779:7568391-7629782:+ 
		ev_e1_s2 = ev_all[2]
		ev_e1_s2 = ev_e1_s2.split('-')
		e1 = ev_e1_s2[0] #discard
		s2 = int(ev_e1_s2[1]) 

		ev_e1_s3 = ev_all[3]
		ev_e1_s3 = ev_e1_s3.split('-')
		s3 = int(ev_e1_s3[1]) 

		if var_region.upper() == "V1":
			s = s3 - len_ext
			motif_coord1 = s + start
			motif_coord2 = s + end
			org_start =  s3 - len_ext
			org_end = s3

		elif var_region.upper() == "V2":
			motif_coord1 = s3 + start
			motif_coord2 = s3 + end
			org_start =  s3
			org_end = s3 + len_ext

		else:
			print("var_region {} not defined for {}. Error!".format(var_region, event))
			quit()

	seq_len = org_end - org_start 
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr,motif_coord1, motif_coord2, 
		event,seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()


def fun_A5_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	##ref: see Fig.2 suppa documentaion of A5 e1,e2 & s3..descriptions
	#https://bitbucket.org/regulatorygenomicsupf/suppa
	#Example event id: (e2-s3:e1-s3)
	#CUL7|9820;A5:chr6:43020534-43021360:43020534-43021449:-
	'''
	fo = open(out_file, "a")
	len_ext = 200
	
	ev_all = event.split(';')	
	ev_coord = ev_all[1]

	ev_all = ev_coord.split(':')
	ev_type = ev_all[0]
	ev_chr = ev_all[1]
	ev_strand = ev_all[4]

	if(ev_strand == "-"):
		s3_e1 = ev_all[2]
		s3_e2 = ev_all[3]
		s3_e1 = s3_e1.split('-')
		s3 = int(s3_e1[0])
		e1 = int(s3_e1[1])

		s3_e2 = s3_e2.split('-')
		e2 = int(s3_e2[1])

		
		if var_region == 'V1':
			s = e2 - len_ext
			motif_coord1 = (e2 - end) 
			motif_coord2 = (e2 - start) 
			org_start =  s
			org_end = e2

		elif var_region == 'V2':
			s = e2 + len_ext
			motif_coord1 = (s - end)
			motif_coord2 = (s - start)
			org_start =  e2
			org_end = s

		else:
			print("var_region not found in event: {}".format(event))

	elif(ev_strand == "+"):
		e2_s3 = ev_all[2]
		e1_s3 = ev_all[3]
		e2_s3 = e2_s3.split('-')
		e2 = int(e2_s3[0])
		s3 = int(e2_s3[1])

		e1_s3 = e1_s3.split('-')
		e1 = int(e1_s3[0])

		if var_region == 'V1':
			s = e1 - len_ext
			motif_coord1 = s + start
			motif_coord2 = s + end
			org_start =  s
			org_end = e1

		elif var_region == 'V2':
			motif_coord1 = e1 + start
			motif_coord2 = e1 + end
			org_start = e1
			org_end = e1 + len_ext
		else:
			print("Error! var_region not found in event: {}".format(event))
			sys.exit()

	else:
		print("No strand information in file, for event: {}".format(event))
		sys.exit()

	#event = event + ";" + motif_id + ";" + kmer
	#gb = ev_chr + ":" + str(motif_coord1) + "-" + str(motif_coord2)

	seq_len = org_end - org_start 
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr,motif_coord1, motif_coord2, 
		event,seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()



def fun_SE_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	ref: see Fig.2 suppa documentaion of SE e1 & s1..descriptions
	https://bitbucket.org/regulatorygenomicsupf/suppa
	Example event id: (e1-s2:e2-s3)
	CTCFL|140690;SE:chr20:56064269-56073610:56073757-56078492:-
	'''
	
	fo = open(out_file, "a")

	#command I am using for fmo conversion to bed
	#path: /projects_rg/babita/TCGA/pipeline/run11/SE/RNAcompete_TCGA/fmo
	#for line in `cat ../motif_homo`; do file=*$line.fmo; ev="SE" ; outf=$ev"."$line".fmo.bed" ; python ~/python_try/fmo_to_bed/test_fmo_to_bed.py $file >../se_fmo_bed/$outf; echo "done: $outf" ; done

	len_ext = 200 #for length extenstion for SE type events  

	ev_all = event.split(';')
	ev_coord = ev_all[1]

	ev_all = ev_coord.split(':')

	try:
		ev_type = ev_all[0]
		ev_chr = ev_all[1]
		ev_e1_s2 = ev_all[2] #intron1_start(e1) to exon_start(s2)
		ev_e2_s3 = ev_all[3] #exon_end(e2) to intron2_end(s3)
		ev_strand = ev_all[4] 

	except Exception, err:
		print("\nError in parsing event_id:\t{}".format(event))
		print("Error:\t{}".format(err))
		sys.exit()


	'''
	double check
	ev_type should be equal to event_type passed from the main function
	'''
	if event_type != ev_type:
		print("Error in event_id: {}, event type {} and {} did not match".format(
				event, event_type, ev_type))
		sys.exit()


	ev_e1_s2 = ev_e1_s2.split('-')

	exon_start_s2 = int(ev_e1_s2[1])  
	ev_e2_s3 = ev_e2_s3.split('-')
	exon_end_e2 = int(ev_e2_s3[0]) 
	intron2_end_s3 = int(ev_e2_s3[1]) 

	if ev_strand == "+":

		if var_region.lower() in ("up","u"):
			s = exon_start_s2 - len_ext
			motif_coord1 = (s + start) - 1
			motif_coord2 = (s + end) - 1
			org_start =  s 
			org_end = exon_start_s2 

		elif var_region.upper() == "E":
			motif_coord1 = (exon_start_s2 + start) -1 #one-base correction previously done
			motif_coord2 = (exon_start_s2 + end) - 1
			org_start = exon_start_s2
			org_end = exon_end_e2

		elif var_region.lower() in ("down", "dow", "do", "d") :
			motif_coord1 = exon_end_e2 + start
			motif_coord2 = exon_end_e2 + end
			org_start = exon_end_e2 
			org_end = exon_end_e2 + len_ext

		else:
			print("var_region {} not defined. Error!".format(var_region))
			sys.exit()
			

	elif ev_strand == '-':

		if var_region.lower() in ("up","u"):
			e = len_ext - start 
			s = len_ext - end
			motif_coord1 = exon_end_e2 + (s+1) 
			motif_coord2 = exon_end_e2 + (e+1)	
			org_start =  exon_end_e2
			org_end = exon_end_e2 + len_ext


		elif var_region.upper() == "E":
			motif_coord1 = (exon_end_e2 - end)  + 1 #one-base correction previously done
			motif_coord2 = (exon_end_e2 - start) + 1
			org_start = exon_start_s2 
			org_end = exon_end_e2 


		elif var_region.lower() in ("down", "dow", "do", "d"):
			motif_coord1 = exon_start_s2 - end
			motif_coord2 = exon_start_s2 - start
			org_start =  exon_start_s2 - len_ext
			org_end = exon_start_s2

		else:
			print("var_region {} not defined. Error!".format(var_region))
			sys.exit()
	else:
		print("No strand information in file, for event: {}".format(event))
		sys.exit()

	seq_len = org_end - org_start 
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr,motif_coord1, motif_coord2, 
		event,seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()


def fun_MX_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	ref: see Fig.2 suppa documentaion of MX e1 & s1..descriptions
	#https://bitbucket.org/regulatorygenomicsupf/suppa
	#Example event id: (e1-s2:e2-s4:e1-s3:e3-s4)
	#SMPD4|55627;MX:chr2:130922018-130925031:130925163-130934127:130922018-130932487:130932573-130934127:-
	#A2BP1|54715;MX:chr16:7703949-7714931:7714970-7726776:7703949-7721559:7721601-7726776:+
	'''
	fo = open(out_file, "a")
	len_ext = 100

	ev_all = event.split(';')	
	ev_coord = ev_all[1]

	ev_all = ev_coord.split(':')
	ev_type = ev_all[0]
	ev_chr = ev_all[1]
	ev_strand = ev_all[6]

	ev_e1_s2 = ev_all[2]
	ev_e2_s4 = ev_all[3]
	ev_e1_s3 = ev_all[4]
	ev_e3_s4 = ev_all[5]

			   
	ev_e1_s2 = ev_e1_s2.split('-') 
	e1 = int(ev_e1_s2[0])
	s2 = int(ev_e1_s2[1])
			  
	ev_e2_s4 = ev_e2_s4.split('-')
	e2 = int(ev_e2_s4[0])
	s4 = int(ev_e2_s4[1])
			  
	ev_e1_s3 = ev_e1_s3.split('-')
	e1 = int(ev_e1_s3[0])
	s3 = int(ev_e1_s3[1])
			  
	ev_e3_s4 = ev_e3_s4.split('-')
	e3 = int(ev_e3_s4[0])
	s4 = int(ev_e3_s4[1])


	#event_mx1 = event + ";MX1"
	#event_mx2 = event + ";MX2"

	if(ev_strand == '+'):

		if var_region == 'V1':
			s = s2 - len_ext
			motif_coord1 = (s + start) - 1
			motif_coord2 = (s + end) - 1
			org_start = s2 - len_ext
			org_end = s2

		elif var_region == 'V2':
			motif_coord1 = (s2 + start) - 1
			motif_coord2 = (s2 + end) - 1
			org_start = s2
			org_end = e2

		elif var_region == 'V3':
			motif_coord1 = e2 + start
			motif_coord2 = e2 + end
			org_start = e2
			org_end = e2 + len_ext

		elif var_region == 'V4':
			s = s3 - 100
			motif_coord1 = (s + start) - 1
			motif_coord2 = (s + end) - 1
			org_start = s3 - len_ext
			org_end = s3

		elif var_region == 'V5':
			motif_coord1 = (s3 + start) - 1
			motif_coord2 = (s3 + end) - 1
			org_start = s3
			org_end = e3

		elif var_region == 'V6':
			motif_coord1 = e3 + start
			motif_coord2 = e3 + end
			org_start = e3
			org_end = e3 + len_ext

		else:
			print("Var_region not found for event:{}".format(event))
			return
			sys.exit()


	elif(ev_strand == '-'):

		if var_region == 'V1':
			s = e3 + len_ext
			motif_coord1 = (s - end) + 1
			motif_coord2 = (s - start) + 1
			org_start = e3 
			org_end = e3 + len_ext

		elif var_region == 'V2':
			motif_coord1 = (e3 - end) + 1
			motif_coord2 = (e3 - start) + 1
			org_start = s3
			org_end = e3

		elif var_region == 'V3':
			motif_coord1 = s3 - end
			motif_coord2 = s3 - start
			org_start = s3 - len_ext
			org_end = s3

		elif var_region == 'V4':
			s = e2 + 100
			motif_coord1 = (s - end) + 1
			motif_coord2 = (s - start) + 1
			org_start = e2
			org_end = s

		elif var_region == 'V5':
			motif_coord1 = (e2 - end) + 1
			motif_coord2 = (e2 - start) + 1
			org_start = s2
			org_end = e2

		elif var_region == 'V6':
			motif_coord1 = s2 - end
			motif_coord2 = s2 - start
			org_start = s2 - len_ext
			org_end = s2

		else:
			print("Var_region not found for event:{}".format(event))
			sys.exit()

	else:
		print("No strand information in file, for event: {}".format(event))
		sys.exit()

	seq_len = org_end - org_start
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr, motif_coord1, motif_coord2, 
		event, seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()


def fun_RI_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	##ref: see Fig.2 suppa documentaion of RI e1 & s1..descriptions
	#https://bitbucket.org/regulatorygenomicsupf/suppa
	#Example event id: 
	#TIAL1|7073;RI:chr10:121336123:121336262-121336592:121336715:-
	'''
	fo = open(out_file, "a")
	len_ext = 200
	
	ev_all = event.split(';')	
	ev_coord = ev_all[1]

	ev_all = ev_coord.split(':')
	ev_type = ev_all[0]
	ev_chr = ev_all[1]
	s1 = int(ev_all[2])
	e1_s2 = ev_all[3]
	e2 = int(ev_all[4])
	ev_strand = ev_all[5]

	e1_s2 = e1_s2.split('-') 
	e1 = int(e1_s2[0])
	s2 = int(e1_s2[1])

	if ev_strand == '+':	
		if var_region == 'V1':
			s = e1 - len_ext
			motif_coord1 = s + start
			motif_coord2 = s + end
			org_start = s
			org_end = e1

		elif var_region == 'V2':
			motif_coord1 = e1 + start
			motif_coord2 = e1 + end
			org_start = e1 
			org_end = s2

		elif var_region == 'V3':
			motif_coord1 = s2 + start
			motif_coord2 = s2 + end
			org_start = s2
			org_end = s2 + len_ext

		else:
			print("No variable region found for event: {}".format(event))
			sys.exit()

	elif ev_strand == '-' :
	
		if var_region == 'V1':
			s = s2 + len_ext
			motif_coord1 = (s - end) 
			motif_coord2 = (s- start)
			org_start = s2
			org_end = s

		elif var_region == 'V2':
			motif_coord1 = s2 - end
			motif_coord2 = s2 - start
			org_start = s2
			org_end = e1

		elif var_region == 'V3':
			motif_coord1 = (e1 - end) + 1
			motif_coord2 = (e1 - start) + 1
			org_start = e1 - len_ext
			org_end = e1

		else:
			print("No variable region found for event: {}".format(event))
			sys.exit()

	else:
		print("No strand information in file, for event: {}".format(event))
		sys.exit()

	#event = event + ";" + motif_id + ";" + kmer
	#gb = ev_chr + ":" + str(motif_coord1) + "-" + str(motif_coord2)

	seq_len = org_end - org_start
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr, motif_coord1, motif_coord2, 
		event, seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()


def fun_AF_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	##ref: see Fig.2 suppa documentaion of A3 e1 & s1..descriptions
	#https://bitbucket.org/regulatorygenomicsupf/suppa
	ABAT|18;AF:chr16:8768444:8768568-8829556:8814573:8814791-8829556:+
	#Example event id: (s1:e1-s3:s2:e2-s3)
	AAA1|404744;AF:chr7:34768428-34797686:34797884:34768428-34800724:34800803:-
	-ve strand: (s3-e2:s2:s3-e1:s1)

	[==V1==]----V2---|--------------------------------[===] AF (+strand)
						[==V3==]----V4-----|----------[===] AF
	'''
	fo = open(out_file, "a")

	len_ext = 200
	ev_all = event.split(';')	
	ev_coord = ev_all[1]

	ev_all = ev_coord.split(':')
	ev_type = ev_all[0]
	ev_chr = ev_all[1]
	ev_strand = ev_all[6]

	if ev_strand == '-':
		#eg.AF:chr7:34768428-34797686:34797884:34768428-34800724:34800803:-
		#(s3-e2:s2:s3-e1:s1)

		ev_s3_e2 = ev_all[2]
		s2 = int(ev_all[3])
		ev_s3_e1 = ev_all[4]
		s1 = int(ev_all[5])

		s3, e2 = ev_s3_e2.split('-')
		s31 = int(s3)
		e2 = int(e2)

		s3, e1 = ev_s3_e1.split('-')
		s32 = int(s3)
		e1 = int(e1)

		if s31 != s32:
			print("Wrong event id. please correct the event: {}".format(event))
			sys.exit()
		
		if var_region == 'V1':
			motif_coord1 = (s1 - end) + 1
			motif_coord2 = (s1 - start) + 1
			org_start = e1
			org_end = s1

		elif var_region == 'V2':
			motif_coord1 = (e1 - end) 
			motif_coord2 = (e1 - start) 
			org_start = e1 - len_ext
			org_end = e1

		elif var_region == 'V3':
			motif_coord1 = (s2 - end ) + 1
			motif_coord2 = (s2 - start) + 1
			org_start = e2
			org_end = s2

		elif var_region == 'V4':
			motif_coord1 = (e2 - end)
			motif_coord2 = (e2 - start)
			org_start = e2 - len_ext 
			org_end = e2

		else:
			print("Error! variable region not found for event: {}".format(event))
			sys.exit()


	elif ev_strand == '+' :  
		#eg. A2BP1|54715;A3:chr16:7568391-7629779:7568391-7629782:+ 
		#(s1:e1-s3:s2:e2-s3)

		s1 = int(ev_all[2])
		ev_e1_s3 = ev_all[3]
		s2 = int(ev_all[4])
		ev_e2_s3 = ev_all[5]

		e1, s3 = ev_e1_s3.split('-')
		e1 = int(e1)
		s31 = int(s3)
		e2, s3 = ev_e2_s3.split('-')
		e2 = int(e2)
		s32 = int(s3)

		if s31 != s32:
			print("Wrong event id. please correct the event: {}".format(event))
			sys.exit()

		if var_region == 'V1':
			motif_coord1 = s1 + start
			motif_coord2 = s1 + end
			org_start = s1 
			org_end = e1

		elif var_region == 'V2':
			motif_coord1 = (e1 + start) + 1
			motif_coord2 = (e1 + end) + 1
			org_start = e1
			org_end = e1 + len_ext

		elif var_region == 'V3':
			motif_coord1 = s2 + start
			motif_coord2 = s2 + end
			org_start = s2 
			org_end = e2

		elif var_region == 'V4':
			motif_coord1 = (e2 + start) - 1
			motif_coord2 = (e2 + end) - 1
			org_start = e2 
			org_end = e2 + len_ext

		else:
			print("Error! variable region not found for event: {}".format(event))
			sys.exit()


	else:
		print("No strand information in file, for event: {}".format(event))
		sys.exit()

	#event = event + ";" + motif_id + ";" + kmer
	#gb = ev_chr + ":" + str(motif_coord1) + "-" + str(motif_coord2)

	seq_len = org_end - org_start
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr, motif_coord1, motif_coord2, 
		event, seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()

def fun_AL_bedfile(event, event_type, var_region, motif_id,
					start, end, kmer, out_file):
	'''
	##ref: see Fig.2 suppa documentaion of A3 e1 & s1..descriptions
	#https://bitbucket.org/regulatorygenomicsupf/suppa
	ABCB5|340273;AL:chr7:20698299-20698689:20700016:20698299-20706382:20707187:+
	#Example event id: (e1-s3:e3:e1-s2:e2)
	ABCA11P|79963:ZNF721|170960;AL:chr4:419224:420762-466364:433781:438221-466364:-
	-ve strand: (e3:s3-e1:e2:s2-e1)

	[===]---------------------------|------V2----[==V1==] AL (+ve strand)
	[===]--------|-----V4-----[==V3==]

	[==V1==]----V2---|--------------------------------[===] AL (-strand)
						[==V3==]----V4-----|----------[===] AL
	'''
	fo = open(out_file, "a")
	len_ext = 200

	ev_all = event.split(';')	
	ev_coord = ev_all[1]

	ev_all = ev_coord.split(':')
	ev_type = ev_all[0]
	ev_chr = ev_all[1]
	ev_strand = ev_all[6]

	if ev_strand == '-':
		#eg.AL:chr4:419224:420762-466364:433781:438221-466364:-
		#(e3:s3-e1:e2:s2-e1)

		e3 = int(ev_all[2])
		ev_s3_e1 = ev_all[3]
		e2 = int(ev_all[4])
		ev_s2_e1 = ev_all[5]

		s3, e1 = ev_s3_e1.split('-')
		s3 = int(s3)
		e11 = int(e1)

		s2, e1 = ev_s2_e1.split('-')
		s2 = int(s2)
		e12 = int(e1)

		if e11 != e12:
			print("Wrong event id. please correct the event: {}".format(event))
			sys.exit()

		if var_region == 'V1':
			motif_coord1 = (s3 - end) + 1
			motif_coord2 = (s3 - start) + 1
			org_start = e3
			org_end = s3

		elif var_region == 'V2':
			s = s3 + len_ext
			motif_coord1 = (s - end) + 1 
			motif_coord2 = (s - start) + 1
			org_start = s3 
			org_end = s

		elif var_region == 'V3':
			motif_coord1 = (s2 - end) + 1
			motif_coord2 = (s2 - start) + 1
			org_start = e2
			org_end = s2

		elif var_region == 'V4':
			s = s2 + 200
			motif_coord1 = (s - end) + 1 
			motif_coord2 = (s - start) + 1
			org_start = s2
			org_end = s

		else:
			print("Error! variable region not found for event: {}".format(event))
			sys.exit()


	elif ev_strand == '+' :  
		#eg. AL:chr7:20698299-20698689:20700016:20698299-20706382:20707187:+
		#(e1-s3:e3:e1-s2:e2)

		ev_e1_s3 = ev_all[2]
		e3 = int(ev_all[3])
		ev_e1_s2 = ev_all[4]
		e2 = int(ev_all[5])

		e1, s3 = ev_e1_s3.split('-')
		e11 = int(e1)
		s3 = int(s3)
		e1, s2 = ev_e1_s2.split('-')
		e12 = int(e1)
		s2 = int(s2)

		if e11 != e12:
			print("Wrong event id. please correct the event: {}".format(event))
			sys.exit()


		if var_region == 'V1':
			motif_coord1 = s3 + start
			motif_coord2 = s3 + end
			org_start = s3 
			org_end = e3

		elif var_region == 'V2':
			s = s3 - 200
			motif_coord1 = (s + start) 
			motif_coord2 = (s + end)
			org_start =  s
			org_end = s3

		elif var_region == 'V3':
			motif_coord1 = s2 + start
			motif_coord2 = s2 + end
			org_start =  s2
			org_end = e2

		elif var_region == 'V4':
			s = s2 - 200
			motif_coord1 = (s + start) 
			motif_coord2 = (s + end)
			org_start =  s
			org_end = s2

		else:
			print("Error! variable region not found for event: {}".format(event))
			sys.exit()


	else:
		print("No strand information in file, for event: {}".format(event))
		sys.exit()


	#event = event + ";" + motif_id + ";" + kmer
	#gb = ev_chr + ":" + str(motif_coord1) + "-" + str(motif_coord2)

	seq_len = org_end - org_start
	fo.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ev_chr, motif_coord1, motif_coord2, 
		event, seq_len,ev_strand, motif_id, kmer, org_start, org_end))
	fo.close()

if __name__ == "__main__":
    main()



