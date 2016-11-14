**Motif Scan and Enrichment Analysis (MoSEA)**
============================================


**MoSEA Modules:**
--------------------
(Please refer test_files/run_test_files.sh to run the pipeline and to run MoSEA directly with SUPPA events).

I. *getfasta*: Extract fasta sequences from given bedfile coordinates.
--------------------------------------------------------------------
	Required Input: 
			 --bedfile : Bed file in atleast 6 columns format, Default:strandness is TRUE
			 --genome : Genome Fasta File sequence in single line, eg. hg19.fa
	Optional Input:
			 --output : Output file name, Default: STDOUT
			 --bedtoolspath : Bedtools bin path, Default: search for 'bedtools' in local path
	
	Output: 1. Sequence Fasta file
	
	Usage: 
	python mosea.py getfasta --bedfile sample.bed --genome /path/to/hg19.fa --output sample_output.fa
	Usage (with optional parameter):
	python mosea.py getfasta --bedfile sample.bed --genome /path/to/hg19.fa --output sample_output.fa --bedtoolspath /path/to/bedtools/bin/bedtools
	
II. *scan*: Scans fasta Sequences for given motifs (PFMs)
-------------------------------------------------
	Required Input: 
			 --scan : Scans for given motif (this option is always true)	

                         --pfm : If selected; Scans for motif using PFMs
                         --pfm_path : Path to folder where PFM matrices of motifs are stored 
                        
                         --kmer : If selected; Scans for motif using k-mers 
                         --kmer_path: Path to k-mer(s) folder, (Note: kmers file should be saved with *.kmer extension, in fasta format or just one kmer in each line. Motif-label will be taken from file_name before ext (eg. for file 'SRSF1.kmer', label will be 'SRSF1').

			 --fasta : Input Fasta file (Sequences in Fasta file format to scan motifs)
			 
			
	Optional Input:
			 --count : Also creates count_file in tabular format for all motifs found on each sequences.
			 --out_dir : Output Directory (path to output directory, Default: creates fmo/ directory in current working path)
			 --fmo_path : Path to FIMO binary folder, Default: Search for 'fimo' in local path
	
	Output: 1. fmo output files for each motif in out_dir/
			2. Summary count_file table for all motifs and all sequences in one file (if --count)}
	
	Usage: 
	python mosea.py scan --pfm --pfm_path /path/to/motif_pfm_dir/ --fasta sample_output.fa
    or
    python mosea.py scan --kmer --kmer_path /path/to/motif_kmer_dir/ --fasta sample_output.fa 

	Usage (with optional parameter):
	python mosea.py scan --scan --pfm --pfm_path /path/to/motif_pfm_dir/  --fasta sample_output.fa  --out_dir /path/my/fmo_output_dir/ --fmo_path /usr/local/bin/fimo --count

	
III. *enrich* : Motif enrichment test analysis between given regulated and background(control) files	
---------------------------------------------------------------------------------------------------
	Required Input: 
			 --reg_fa_file : Input file for regulated sequences in fasta format
			 --reg_count_file : Input file for regulated motif counts in tab seperated format(header=Motif_ids)
			 --bg_fa_file : Input file for background (control) sequences in fasta format
			 --bg_count_file : Input file for background (control) motif counts in tab seperated format(header=Motif_ids)
			
	Optional Input:
			 --out_file : Creates output for enrichment test analysis, Default: output_MoSEA_enrichment.txt
			 --rand_count : Number of times randomization performed for each regulated sequence. Default:100
			 --match_len : Match length of background sequences to each regulated sequence(1=True, 0=False). Default: 0 (False)
			 --len_ext : Percent range of length extension on both side boundaries, to match with background sequences.
			 			 This option is valid only when '--match_len 1'. 
						 Default: 15 (15% length extended on both sides)
	
	Output: 1. Enrichment test output file for each motif with zscore 
	
	Usage: 
	python mosea.py enrich --reg_fa_file reg_sample_output.fa --reg_count_file reg_sample_output.space
						   --bg_fa_file bg_sample_output.fa --bg_count_file bg_sample_output.space
						 
	Usage (with optional parameter):
	python mosea.py enrich --reg_fa_file reg_sample_output.fa --reg_count_file reg_sample_output.space
						   --bg_fa_file bg_sample_output.fa --bg_count_file bg_sample_output.space
						   --out_file my_output_file.tab --rand_count 200 --match_len 1 --len_ext 20
						   



_Run MoSEA on Suppa events to get enriched Motifs_
-------------------------------------------------
Requirements:
```
- 	Python 2.7 (numpy, bio, pandas, StringIO, collections, math)
- 	Bedtools
- 	Fimo (MEME-Suite version >version 4.11 or higher)
- 	Genome fasta sequence (\*.fa and \*.fai files)
- 	Motif PFMs or Kmers file or directory 
```	
Input files:
````
- 	Suppa event ids for regulated set (or Bedfile cordinates)
- 	Suppa event ids for control set (or Bedfile cordinates)
```	

The analysis is divided into two steps, A) Motif Scan B) Enrichment

_For Analysis (A): Motif Scan Steps:_
-------------------------------------------------
Input: 2 files: Regulated & Control event ids (Suppa events or Bedfile format)
Steps: Input_file -> Extract Sequence -> Scan sequences for Motif -> Create count table
Requirement: Genome Fasta file (hg19.fa), Bedtools, FIMO from MEME suite

_For Analysis (B): Enrichment_
-------------------------------------------------
Input: 4 files: (All these files are generated in Analysis A)
- 	Regulated fasta sequences file 
-  	Regulated motif count table (generated in step A)
-  	Control(Background) fasta sequence file
-  	Control motif count file (generated in step A)

Steps: For each regulated sequence -> create pools of Control seqs matching GC content and length of seq -> randomize 100 times -> count motifs on reg & control -> calculate z-score by observe(reg) vs expected (distribution from control) ((x-mean)/SD) 

To test run call this script from MoSEA directory: 
cd MoSEA/
./test_files/run_test_files.sh (Please point the paths to appropriate directories)


Commands:
-----------------------------------------------------------------------------------------------------

## Analysis (A): Motif Scan 


#### Step1 : Convert Suppa events to bedfile format (This step is not necessary if cordinates are already in bedfile format. Go to step 2 directly.) 
```
$path_python ./mosealib/suppa_to_bed.py --ifile $reg_infile --event SE --ext 200 --ofile $bedfile_reg
```
(Ext : Extension is upstream and downstream regions from the center Exon)

![Extension img](https://github.com/comprna/MoSEA/blob/master/img/var_regions_suppa.jpg)


#### Step2:  Convert Bedfile to Fasta File
```
$path_python mosea.py getfasta --bedfile $bedfile_reg --genome $path_genome --output $fafile_reg
```

#### Step3:   Scan fasta sequences for Motifs (Example shown for PFMs)
```
fmopfm_outdir="fmo_pfm"  #output dir for scanned motifs
$path_python mosea.py scan --pfm --pfm_path $path_pfms --fasta $fafile_reg --out_dir $fmopfm_outdir --fmo_path $path_fimo --count
```

Repeat above three steps for control files as well. Ideally, number of events in control set must be atleast 100x more than regulated set.


## Analysis (B): Motif Enrichment Analysis 


Required Input:  4 files (2 regulated files: fasta & count_table, 2 control files: fasta & count_table) generated in Analysis (A).

#### Step4: Perform Enrichment
```
Regulated files : Fasta sequences and Motif count table
- reg_file_fa
- reg_file_count
	
Control files : Fasta sequences and Motif count table
- control_file_fa
- control_file_count

#Zscore Output file name
outfile="zscore_outfile.tab"

#perform enrichment
$path_python mosea.py enrich --reg_fa_file $reg_file_fa --reg_count_file $reg_file_count \
                       --bg_fa_file $control_file_fa --bg_count_file $control_file_count \
		     	--out_file $outfile

```
-------------------------------------------------
Please run the dummy script in test_files/run_test_files.sh to view the outputs at each step.
