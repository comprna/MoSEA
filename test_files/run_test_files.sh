event="SE"
eventl=${event,,} #lower case
out_dir="/projects_rg/babita/TCGA/soft/mosea/test_files"
hg19_genome="/projects_rg/babita/TCGA/miRNA/db/genome/hg19.fa"
mosea="/projects_rg/babita/TCGA/soft/mosea"
dpsi_file="/projects_rg/TCGA/pipeline/run11/brca_event_dpsi-"$eventl"_paired-filtered.txt"
pfm_path="/home/babita/tmp/pfm_tmp"
python="/soft/devel/python-2.7/bin/python"

#create events
cat $dpsi_file | sed '1d' | cut -f2 | head -80 >"$out_dir/test_"$event"_events"
#convert events to bedfiles
$python $mosea/mosealib/create_bedfiles_from_suppa_events.py --ifile "$out_dir/test_"$event"_events" --event $event --ext 200 --ofile "$out_dir/test_"$event"_events.bed"
#convert bedfiles to fasta files
$python $mosea/mosea.py getfasta --bedfile "$out_dir/test_"$event"_events.bed" --genome $hg19_genome --output "$out_dir/test_"$event"_events.fa"
#scan motifs on fasta file
$python $mosea/mosea.py scan --scan --count --fasta "$out_dir/test_"$event"_events.fa" --pfm_path $pfm_path --out_dir $out_dir
#renmae fmo file
cat "$out_dir/test_"$event"_events.fa.motif_00001.fmo" >"$out_dir/test_"$event"_motif_00001.fmo"
rm "$out_dir/test_"$event"_events.fa.motif_00001.fmo"




