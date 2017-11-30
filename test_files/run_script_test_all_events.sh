#Pipeline to run MoSEA (Motif enrichment Analysis)
#Input: Suppa event ids or cordinates in Bedfile format.

#Note: The analysis is divided into two steps, A) Motif Scan B) Enrichment
#For Analysis (A): Motif Scan Steps:
#====================================
#Input: 2 files: Regulated & Control event ids (Suppa events or Bedfile format)
#Steps: Input_file -> Extract Sequence -> Scan sequences for Motif -> Create count table
#Requirement: Genome Fasta file (hg19.fa), Bedtools, FIMO from MEME suite

#For Analysis (B): Enrichment
#====================================
#Input: 4 files: (All these files are generated in Analysis A)
#- Regulated fasta sequences file 
#- Regulated motif count table 
#- Control(Background) fasta sequence file
#- Control motif count file

#Steps: For each regulated sequence -> create pools of Control seqs matching GC content and length of seq -> randomize 100 times -> count motifs on reg & control -> calculate z-score by observe(reg) vs expected (distribution from control) ((x-mean)/SD) 

#Call this script after defining following paths: 
#./run_script_test_all_events.sh 

#==============================================================
# Important:
#==============================================================
# *** This script is for extended analysis for all Suppa events, for simple generic analysis please see script 'run_test_files.sh' in mosea/test_files directory. ***
# Search '#change' in the text to change the paths required. (rest is optional, the script will generate the filenames itself by provided prefix).


#change pathfor output dir
path_outdir="./mosea_files" 			#change

#define paths:
path_python="/soft/devel/python-2.7/bin/python"	#change
path_mosea="/path/to/mosea"			#change
path_genome="/path_to/genome/hg19.fa"		#change

#not need to change following paths, although please have a look.
fmo_outdir="$path_outdir/fmo"
fmobed_outdir="$path_outdir/fmobed"
path_pfms="$path_mosea/test_files/motifs/pfms"
path_kmers="$path_mosea/test_files/motifs/kmer"
path_fimo="$path_mosea/test_files/motifs/fimo"

mkdir $path_outdir
mkdir $fmo_outdir
mkdir $fmobed_outdir


create_events_and_scan_motif() {

se_upper=$1
events_infile=$2
file_prefix=$3 #prefix of outfiles
		   
 if [ $se_upper == "MX" ] #for MX events we usually take 100ntd variable regions
 then
    len_ext=100
 else
    len_ext=200
 fi 
 
 for reg in "${regulation[@]}" ; do

     se_lower=${se_upper,,}  #convert to lower case

     outfile=$file_prefix"_"$se_lower"_"$reg	 #outfile name
     bedfile=$path_outdir"/"$outfile".bed"	
     fastafile=$path_outdir"/"$outfile".fa"

     #convert events to bedfile
     $path_python $path_mosea/mosealib/create_bedfiles_from_suppa_events.py \
	    	--ifile "$events_infile" --event $se_upper --ext $len_ext \
	    	--ofile $bedfile

     #convert bedfile to fasta
     $path_python $path_mosea/mosea.py getfasta --bedfile $bedfile \
	    	--genome $path_genome --output $fastafile

      #scan fasta sequences
      $path_python $path_mosea/mosea.py scan --scan --count --fasta $fastafile \
	     --out_dir $fmo_outdir --pfm_path $path_pfms --fmo_path $path_fimo


 done
}

perform_motif_enrichment(){
#Call python script for enrichment analysis using MoSEA

	$path_python $path_mosea/mosea.py enrich --reg_fa_file $1 \
		--reg_count_file $2 --bg_fa_file $3 --bg_count_file $4 \
		--out_file $5

}

create_files_for_enrichment() {

se_upper=$1
file_prefix=$2

mkdir tmp

  for reg in "${regulation2[@]}" ; do
	   
	if [ $se_upper = "SE" ]
	then  			
	location=("up" "E" "down")
		   
	elif [ $se_upper == "MX" ]
	then
	location=("V1" "V2" "V3" "V4" "V5" "V6") 
	  		
	elif [ $se_upper == "RI" ]
	then
	location=("V1" "V2" "V3") 

	elif [ $se_upper == "A3" ]
	then
	location=("V1" "V2") 

	elif [ $se_upper == "A5" ]
	then
	location=("V1" "V2") 

	fi 
		
	 #call SE function with variable region
	#		 for loc in "${location[@]}" ; do
	
	for loc in "${location[@]}" ; do 
			
	se=${se_upper,,}
	reg_file_fa=$file_prefix"_"$se"_"$reg".fa"
	reg_file_count=$file_prefix"_"$se"_"$reg".fa.tab"	
	control_file_fa=$file_prefix"_"$se"_control.fa"
	control_file_count=$file_prefix"_"$se"_control.fa.tab"
	outfile=$file_prefix"_"$se"_"$reg"-"$loc"_motif_zscore.tab"
	pat=";"$loc
				
	echo "For: $se | $file_prefix | $reg | $loc "
	grep -A1 "$pat" $path_outdir/$reg_file_fa >tmp/tmp_reg_file_fa
	grep -A1 "$pat" $path_outdir/$control_file_fa >tmp/tmp_cont_file_fa

	cat $fmo_outdir/$reg_file_count | head -1 >tmp/tmp_reg_file_count; grep "$pat" $fmo_outdir/$reg_file_count >>tmp/tmp_reg_file_count
	cat $fmo_outdir/$control_file_count | head -1 >tmp/tmp_cont_file_count; grep "$pat" $fmo_outdir/$control_file_count >>tmp/tmp_cont_file_count
				
	perform_motif_enrichment tmp/tmp_reg_file_fa tmp/tmp_reg_file_count tmp/tmp_cont_file_fa tmp/tmp_cont_file_count $outfile
	
	rm tmp/tmp_reg_file_fa tmp/tmp_reg_file_count tmp/tmp_cont_file_fa tmp/tmp_cont_file_count
   		
	done	
	 
   done

rm -rf tmp
}


#Run script for MoSEA enrichment analysis

regulation=("control" "negative" "positive")
regulation2=("positive" "negative")
se_type=("RI" "A3" "MX" "A5" "SE" )


#create bed file.
#input file, se_type, outdir

##call functions
for se in "${se_type[@]}" ; do
	
	 input_event_file="$se"_list_of_events.txt 	#change 
	 file_prefix="file" 				#change  (Added before each file name)

	 echo "start for $se "
	 create_events_and_scan_motif $se $input_event_file $file_prefix
	 create_files_for_enrichment $se $file_prefix

done



	


