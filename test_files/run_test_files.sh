#Pipeline to run MoSEA (Motif enrichment Analysis)
#Input: Suppa event ids or cordinates in Bedfile format.

#Note: The analysis is divided into two steps, A) Motif Scan B) Enrichment
#For Analysis (A): Motif Scan Steps:
#====================================
#Input: One File with Suppa events or Bedfile format  
#Steps: Input_file -> Extract Sequence -> Scan sequences for Motif -> Create count table
#Requirement: Genome Fasta file (hg19.fa), Bedtools, FIMO from MEME suite

#For Analysis (B): Enrichment
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

path_fimo="~/meme/bin/fimo"
path_genome="./MoSEA/test_files/genome/hg19_chr22.fa"
path_mosealib="./MoSEA/mosealib"
path_pfms="./MoSEA/test_files/pfms"
path_kmers="/MoSEA/test_files/kmer"
path_outdir="./MoSEA/test_files/outdir"



#Analysis (A): Motif Scan
#====================================

#-------------------------------------------------------
#1. For Regulated Suppa events
#-------------------------------------------------------

infile="./MoSEA/test_files/infile/reg_events_chr22.ids"
outfile_reg="$path_outdir/reg_events_chr22.bed"

#Step1: Convert Suppa events to bedfile format (This step is not necessary if cordinates are already in bedfile format. Go to step 2 directly.) 
python $path_mosealib/suppa_to_bed.py --ifile "$infile" --event SE --ext 200 --ofile $outfile_reg


#Step2: Convert Bedfile to Fasta File
infile="$path_outdir/$outfile1" #output file created in Step1:
reg_outfile_fa="$path_outdir/$infile.fa"
python mosea.py getfasta --bedfile $outfile --genome $path_genome --output $reg_outfile_fa

   
#Step3:  #Scan fasta sequences for Motifs
infile="$path_outdir/$outfile_fa" #output file created in Step2:
fmopfm_outdir="$path_outdir/fmo_pfm"   

#Scan for PFMs
python mosea.py scan --pfm --pfm_path $path_pfms --fasta $reg_outfile_fa --out_dir $fmopfm_outdir --fmo_path $path_fimo --count


#-------------------------------------------------------
#Repeat the above steps 1,2,3 for Control file:
#--------------------------------------------------------

#For Control Suppa events
infile="./MoSEA/test_files/infile/control_events_chr22.ids"
outfile_control="$path_outdir/control_events_chr22.bed"
python $path_mosealib/suppa_to_bed.py --ifile "$infile" --event SE --ext 200 --ofile $outfile_control

control_outfile_fa="$path_outdir/$outfile1.fa"
python mosea.py getfasta --bedfile $outfile1 --genome $path_genome --output $control_outfile_fa

#Scan fasata sequqnces for PFMs
python mosea.py scan --pfm --pfm_path $path_pfms --fasta $control_outfile_fa --out_dir $path_outdir/fmo_pfm --fmo_path $path_fimo --count




#========================================
#Analysis (B): Motif Enrichment Analysis
#========================================

#Required Input:  4 files (2 regulated files: fasta & count_table, 2 control files: fasta & count_table) generated in Analysis (A).


#Regulated files : Fasta sequences and Motif count table
reg_file_fa=$reg_outfile_fa
reg_file_count="$reg_file_fa.tab"	
	
#Control files : Fasta sequences and Motif count table
control_file_fa="$control_outfile_fa"    
control_file_count="$control_file_fa.tab"

#Zscore Output file name
outfile="$path_outdir/zscore_outfile.tab"

#perform enrichment
python mosea.py enrich --reg_fa_file $reg_file_fa --reg_count_file $reg_file_count \
                       --bg_fa_file $control_file_fa --bg_count_file $control_file_count \
		       --out_file $outfile

	






	




