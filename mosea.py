import argparse
import sys,os
#custom libraries
from mosealib import controller
from mosealib import enrichment

#Things to do:
#1. Include pfm matrix creator
def main():
    #create a top-level parser
    parser = argparse.ArgumentParser(description = "Runs motif scan & enrichment")
    subparsers = parser.add_subparsers(help="commands help")
    
    #create parser for 'Get Fasta'
    #parser_fasta = argparse.ArgumentParser(description = "Extract Fasta file from Bed file")
    #parser_fasta = subparsers.add_parser(help="comamnds")
    fasta_parser = subparsers.add_parser("getfasta", 
        help="Extract Fasta sequences from input bedfile and genome.fa")
    fasta_parser.add_argument("--bedfile", "-b", 
        help="Coordinates to  extract in bed format")
    fasta_parser.add_argument("--genome", "-g", metavar="GENOME_FASTA_FILE",
        help="Single Genome fasta file")
    fasta_parser.add_argument("--output", "-o", 
        default=sys.stdout, help="output fasta file name")
    fasta_parser.add_argument("--bedtoolspath", "-path", default="bedtools", 
        help="Bedtools location (if its not in default path)")
    fasta_parser.set_defaults(func=getfasta)


    # Scan Sequence by PFM & Kmers
    scan_parser = subparsers.add_parser("scan", 
        help="Perform motif scanning on fasta sequences")
    scan_parser.add_argument("--pfm", default=False, 
        action="store_true", help="Scan Motif PFMs file")
    scan_parser.add_argument("--pfm_path",
        help="Path to PFMs folder")
    scan_parser.add_argument("--fmo_path","-fmopath", default="fimo",
        help="Path to FIMO binary folder (if not in default path)")
    scan_parser.add_argument("--kmer", default=False, 
        action="store_true", help="Scan Motif Kmers file")
    scan_parser.add_argument("--kmer_path", 
        help="Path to Kmers folder")

    scan_parser.add_argument("--count", default=False, 
        action="store_true", help="create count-Motif file")
    scan_parser.add_argument("--fasta", "-fa", 
        help="Input Fasta file name to scan motifs")
    scan_parser.add_argument("--out_dir","-odir", default="fmo",
        help="Output Folder")

    scan_parser.set_defaults(func=scan)



    #enrichment
    enrich_parser = subparsers.add_parser("enrich",
        help="Perform enrichment analysis between regulated and control\
        (background) set")
    enrich_parser.add_argument("--reg_fa_file", default=False,
         help="input file for regulated sequences in\
        fasta format")
    enrich_parser.add_argument("--reg_count_file", default=False,
         help="input file for regulated motif counts\
        file in table format(header=Motif_ids)")
    enrich_parser.add_argument("--bg_fa_file", default=False,
         help="input file for background(control) \
        sequences in fasta format")
    enrich_parser.add_argument("--bg_count_file", default=False,
         help="input file for background(control) \
        motif counts file in table format(header=Motif_ids)")
    enrich_parser.add_argument("--out_file", default="output_MoSEA_enrichment.txt",
        help="output filename for enrichment analysis")
    enrich_parser.add_argument("--rand_count", metavar='int', default=100,
        type=int, help="Number of times randomization performed \
        for each sequence. Default:100")
    enrich_parser.set_defaults(func=enrich)


    # Overlap
    overlap_parser = subparsers.add_parser(
        "overlap", help="Compares motif lists.")
    overlap_parser.add_argument("input_file")
    overlap_parser.add_argument(
        "--strict", "-s",  default=False, action="store_true")
    overlap_parser.set_defaults(func=overlap)
    

    try:
        args = parser.parse_args()
        args.func(args)

    # except Exception, err:
    #     print("\nError in parsing: {}".format(err))
    #     print(parser.format_help())
    #     sys.exit(1)

    except Exception as error:
        print("\nError in parsing: {}".format(error))
        print(parser.format_help())
        sys.exit(1)

def getfasta(args):

    controller.get_fasta_sequence(args.bedtoolspath, args.bedfile, args.genome, 
            args.output)
    print("Extracted Fasta sequence saved in: {}".format(args.output))
    
def scan(args):

    if args.pfm:
        controller.scan_motif_pfm(args.fasta, args.out_dir, args.pfm_path, args.fmo_path, args.count)

    if args.kmer:
        controller.scan_motif_kmer(args.fasta, args.out_dir, args.kmer_path, args.count)
        
    #if args.count:
     #   controller.count_motif(args.fasta, args.pfm_path, args.out_dir)  
          
def enrich(args):
    enrichment.calculate_enrichment(args.reg_fa_file, args.reg_count_file,
                                    args.bg_fa_file, args.bg_count_file,
                                    args.out_file, args.rand_count)

def overlap(args):
    print("Looking for overlaps!")
    #print(args.strict)



    
if __name__ == "__main__":
    main()



#pwd: /home/babita/Documents/scan_motifs
#command 1: python scripts/find_motif_overlaps.py getFasta -b mock_files/tmp_bed_file -g mock_files/hg19_genome.fa -o output/tmp_bed_file.fa
#command 2: python scripts/find_motif_overlaps.py scan -fa output/tmp_bed_file.fa -odir output/fmo/ -pfm mock_files/pfm/

