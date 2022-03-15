### no-thinking workflow for processing bining data ###
### code author: Yufei Zeng
### yfzeng0827@hotmail.com
### github.com/Yflyer/

############### before the work ##################
# please make sure complete the previous tutorial 
# enter the project directory
mkdir Bin
cd Bin

# make a work screen
# screen is useful, please know more
screen -s binning

# We still use one sample to complete whole workflow.
# If you want to process multi samples, please refer to parallel tools or loop programing
ln -s ../02_megahit/*/*.fa .
ln -s ../01_bbmap/*clean.fq .

##################  Preparation  ################
### individual assembly stratgy
conda activate py36
### most of binning tools need contigs > 2000 bps 
### if your samples are too fragmented, you can set its length lower (e.g., 1000).)
seqkit seq -m 1500 S1.contigs.fa > S1.m1500.fa

### Binning needs abundance data to distinguish reads from the same host
### we used coverm to retrive contigs abundance
conda activate coverM
coverm contig --interleaved *.clean.fq --reference S1.m1500.fa -t 60 --bam-file-cache-directory contig_bam -o S1.cov.tsv

##################  Binner use  ################
### Tool: maxbin
conda activate py36
mkdir S1.maxbin
### We need to modify abundance file to fit in maxbin format. Please refer more details to maxbin's README
# file excluding first N lines: tail -n +<N+1> <filename>
cat S1.cov.tsv | tail -n +2 | awk '{print $1"\t"$2}' > S1.maxbin.cov.tsv
run_MaxBin.pl -contig S1.m1500.fa -abund S1.maxbin.cov.tsv -out S1.maxbin/S1 -thread 6

### Tool: metabat
conda activate metabat2
jgi_summarize_bam_contig_depths --outputDepth S1.metabat.cov.tsv contig_bam/S1*.bam
metabat2 -i S1.m1500.fa -a depth.txt -o S1.metabat/S1

### Now we finished the use of binner, then we used das tool to filter
### (OPTIONAL:)Or you can use metawrap pipeline to run it once
conda acitvate metawrap
metawrap binning -o ${i}.metawrap -t 20 -a ${i}/${i}.contigs.fa --metabat2 --maxbin2 --concoct --interleaved ${i}/*.fastq
###  Unaddressed method: Vamb

##################  Quality control  ################
### das_genome to auto-imporve genome's quality
conda activate 
# generate dastool-used contig list of each method
Fasta_to_Scaffolds2Bin.sh -i S1.maxbin/ -e fasta > S1.maxbin2das.tsv
Fasta_to_Scaffolds2Bin.sh -i S1.metabat/ -e fa > S1.metabat2das.tsv
# input all the result by each method
mkdir S1.dastool
DAS_Tool -i S1.metabat2das.tsv,S1.maxbin2das.tsv --write_bins 1 -l metabat,maxbin -c S1.m1500.fa -o S1.dastool/S1 -t 15

### check M
# pplacer_threads: number of threads used by pplacer (memory usage increases)
# nt: generate nucleotide gene sequences for each bin
# (out).0XX.fasta	the XX bin. XX are numbers, e.g. out.001.fasta
# (out).summary	summary file describing which contigs are being classified into which bin.
# (out).log	log file recording the core steps of MaxBin algorithm
# (out).marker	marker gene presence numbers for each bin. This table is ready to be plotted by R or other 3rd-party software.
# (out).marker.pdf	visualization of the marker gene presence numbers using R
# (out).noclass	all sequences that pass the minimum length threshold but are not classified successfully.
# (out).tooshort	all sequences that do not meet the minimum length threshold.

### lineage workflowauto-run tree, lineage_set, analyze, qa
checkm lineage_wf --nt -x fa -t 8 --pplacer_threads 4 --tab_table -f S1.checkm.tsv S1.dastool/S1_DASTool_bins S1.checkm

### (OPTIONAL:)check M for CPR:
# Identify marker genes in bins and calculate genome statistics
checkm analyze -x fa -t 8 ~/hmm_sets/cpr_43_markers.hmm S1.dastool/S1_DASTool_bins S1.checkm.cpr
# Assess bins for contamination and completeness.
checkm qa -t 8 -o 1 -f test.cpr.tsv --tab_table ~/hmm_sets/cpr_43_markers.hmm S1.checkm.cpr
