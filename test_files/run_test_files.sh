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

#Call this script from MoSEA directory: 
#cd MoSEA/
#./test_files/run_test_files.sh 


#Requirements: 
#Change following paths if necessary

path_python="/soft/devel/python-2.7/bin/python" #Script not checked for Python 3
path_fimo="./mosealib/fimo" #Otherwise fimo from MEME-Suite (version 4.11 or higher, Careful:older version truncates Suppa event_ids)
path_genome="./test_files/genome/hg19.fa" #file should have index file *.fai in same directory (eg: hg19.fa.fai)
path_pfms="./test_files/motifs/pfms"	  #It could be a dir or just one *.pfm file (ALL pfm files should be with *.pfm extention)
path_kmers="./test_files/motifs/kmer"	#could be a path to dir or path to one file; all files should be with *.kmer extension

#Change Input file parameters here: regulated events, control events
#---------------------------------------------------------
reg_infile="./test_files/dummy/infile/reg_events.ids"
control_infile="./test_files/dummy/infile/control_events.ids"
path_outdir="./test_files/dummy/outfile"
#---------------------------------------------------------



#Analysis (A): Motif Scan
#====================================

#-------------------------------------------------------
#1. For Regulated Suppa events
#-------------------------------------------------------

fname=$(basename $reg_infile)
bedfile_reg="$path_outdir/$fname.bed"

#Step1: Convert Suppa events to bedfile format (This step is not necessary if cordinates are already in bedfile format. Go to step 2 directly.) 
$path_python ./mosealib/suppa_to_bed.py --ifile $reg_infile --event SE --ext 200 --ofile $bedfile_reg


#Step2: Convert Bedfile to Fasta File
fafile_reg="$bedfile_reg.fa"  #output file name
$path_python mosea.py getfasta --bedfile $bedfile_reg --genome $path_genome --output $fafile_reg
 
#Step3:  #Scan fasta sequences for Motifs (Example shown for PFMs)
fmopfm_outdir="$path_outdir/fmo_pfm"  #output dir for scanned motifs
$path_python mosea.py scan --pfm --pfm_path $path_pfms --fasta $fafile_reg --out_dir $fmopfm_outdir --fmo_path $path_fimo --count

#-------------------------------------------------------
# For Control Suppa events (Repeat above steps 1,2,3):
#--------------------------------------------------------
fname=$(basename $control_infile)
bedfile_control="$path_outdir/$fname.bed"


#Covert ids to bedfile
$path_python ./mosealib/suppa_to_bed.py --ifile $control_infile --event SE --ext 200 --ofile $bedfile_control

#Extract fasta sequence
fafile_control="$bedfile_control.fa"
$path_python mosea.py getfasta --bedfile $bedfile_control --genome $path_genome --output $fafile_control

#Scan fasta sequences for PFMs
$path_python mosea.py scan --pfm --pfm_path $path_pfms --fasta $fafile_control --out_dir $path_outdir/fmo_pfm --fmo_path $path_fimo --count




#========================================
#Analysis (B): Motif Enrichment Analysis
#========================================

#Required Input:  4 files (2 regulated files: fasta & count_table, 2 control files: fasta & count_table) generated in Analysis (A).


#Regulated files : Fasta sequences and Motif count table
reg_file_fa=$fafile_reg
reg_file_count="$fafile_reg.tab"
	
#Control files : Fasta sequences and Motif count table
control_file_fa=$fafile_control  
control_file_count="$fafile_control.tab"

#Zscore Output file name
outfile="$path_outdir/zscore_outfile.tab"

#perform enrichment
$path_python mosea.py enrich --reg_fa_file $reg_file_fa --reg_count_file $reg_file_count \
                       --bg_fa_file $control_file_fa --bg_count_file $control_file_count \
		       --out_file $outfile

	






	




