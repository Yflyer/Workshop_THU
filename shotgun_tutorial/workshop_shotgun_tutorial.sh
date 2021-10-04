### no-thinking workflow for processing single metagenomic data ###
### code author: Yufei Zeng
### yfzeng0827@hotmail.com
### github.com/Yflyer/

############### before the work ##################
# please copy the data to your project directory
mkdir [your_project_name]
cd [your_project_name]
cp -r /vd04/yufei/workshop_shotgun/test/0_rawdata [your_project_name]
# make a work screen
# screen is useful, please know more
screen -s workshop

############ steps in our workshop ###############
### The goal of workflow is to generate a GeoChip-like result ###
# quality control (qc)
# contigs assembly
# ORFs prediciton
# sequence mapping
# ORFs cluster (we will get a GeoChip-like result at this step)
# bining
# (annotation, uncompleted)

#### activate conda env #####
conda activate py36
###################  qc  ######################
# time_consuming: ★☆
mkdir 01_cleandata
cd 01_cleandata
ln -s ../0_rawdata/* ./ # link
trimmomatic PE -phred33 -threads 4 S1_r1.fq.gz S1_r2.fq.gz \
      trimmed.S1_r1.fq.gz outtrimmed.S1_r1.fq.gz trimmed.S1_r2.fq.gz outtrimmed.S1_r2.fq.gz  \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 \
      MINLEN:50

##################   contig assembly   #######################
### we used megahit to save time and memory
# time_consuming: ★★★★
cd ..
mkdir -p 02_megahit
cd 02_megahit
ln -s ../01_cleandata/trimmed* ./
megahit -1 trimmed.S1_r1.fq.gz -2 trimmed.S1_r2.fq.gz --min-count 2 --k-list 29,39,51,67,85,107,133 -m 0.5 -t 4 --min-contig-len 500 --out-prefix S1 -o S1

sed -i "s/>/>S1\_/1" S1/S1.contigs.fa

#################   ORF prediciton  #######################
### we used prokka to avoid too much detail
# time_consuming: ★★☆
cd ..
mkdir -p 03_prokka
cd 03_prokka
ln -s ../02_megahit/*/*.fa ./
prokka --metagenome --cpus 4 --outdir S1 --prefix S1 --mincontiglen 500 S1.contigs.fa

sed -i "s/>/>S1\_/1" S1/S1.ffn
sed -i "s/>/>S1\_/1" S1/S1.faa
############## prokka
#.gff	This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV.
#.gbk	This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence.
#.fna	Nucleotide FASTA file of the input contig sequences.
#.faa	Protein FASTA file of the translated CDS sequences.
#.ffn	Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
#.sqn	An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc.
#.fsa	Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines.
#.tbl	Feature Table file, used by "tbl2asn" to create the .sqn file.
#.err	Unacceptable annotations - the NCBI discrepancy report.
#.log	Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled.
#.txt	Statistics relating to the annotated features found.
#.tsv	Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product

##################   mapping  ################
### we used bamm to quickly get mapping and coverage
# time_consuming: ★★★☆
cd ..
#switch to conda env py27
conda activate py27
mkdir -p 04_mapping
cd 04_mapping
ln -s ../03_prokka/*/*.ffn ./
ln -s ../01_cleandata/trimmed* ./

bamm make -d S1.ffn -c trimmed.S1_r1.fq.gz trimmed.S1_r2.fq.gz -t 4 --out_folder S1
bamm parse -c S1.covs.tsv -b S1/S1.ffn.trimmed.S1_r1.bam

#### Coverage calculation modes
### BamM implements several coverage calculation methods. The user can choose the
### method using the -m argument. ###
# opmean: Outlier pileup coverage: average of reads overlapping each base, after bases
# with coverage outside mean +/- 1 standard deviation have been excluded. The number of
# standard deviation used for the cutoff can be changed with --coverage_range.
# pmean: Pileup coverage: average of number of reads overlapping each base
# tpmean: Trimmed pileup coverage: average of reads overlapping each base, after bases
# with in the top and bottom 10% have been excluded. The 10% range can be changed
# using --coverage_range.
# counts: Absolute number of reads mapping
# cmean: Like 'counts' except divided by the length of the contig
# pmedian: Median pileup coverage: median of number of reads overlapping each base

##################  cluster   ################
### we used mmseq to quickly cluster
# time_consuming: ★
cd ..
conda activate py36
mkdir 05_orf_clustering
cd 05_orf_clustering

cat ../02_megahit/*/*.fa >> merge.contigs.fa
cat ../03_prokka/*/*.faa >> merge.faa
cat ../03_prokka/*/*.ffn >> merge.ffn
# cluster and get its info
mkdir db clu_db clu_rep

# contigs cluster
mmseqs createdb merge.contigs.fa db/contigs
mmseqs linclust db/contigs clu_db/contigs tmp --threads 4 --min-seq-id 0.8
mmseqs createtsv db/contigs db/contigs clu_db/contigs clu.contigs.tsv --threads 4
mmseqs createseqfiledb db/contigs clu_db/contigs clu_rep/contigs
mmseqs result2flat db/contigs db/contigs clu_rep/contigs clu_rep.contigs.fasta # get fasta file

# faa cluster
mmseqs createdb merge.faa db/faa
mmseqs linclust db/faa clu_db/faa tmp --threads 4 --min-seq-id 0.8
mmseqs createtsv db/faa db/faa clu_db/faa clu.faa.tsv --threads 4
mmseqs createseqfiledb db/faa clu_db/faa clu_rep/faa
mmseqs result2flat db/faa db/faa clu_rep/faa clu_rep.faa.fasta # get fasta file
#################################################

### and next, we change to R environment to process our ORFs data

### Bonus: bin your contigs, to get metegenomic-assembly genomes (MAGs) ###
# time_consuming: ★★★★
cd ..
mkdir 06_bining
cd 06_bining
### first: get contigs coverage as abundunce file
ln -s ../02_megahit/*/*.fa ./
ln -s ../01_cleandata/trimmed* ./

conda activate py27
bamm make -d S1.contigs.fa -c trimmed.S1_r1.fq.gz trimmed.S1_r2.fq.gz -t 4 --out_folder S1
bamm parse -c S1.covs.tsv -b S1/S1.contigs.trimmed.S1_r1.bam

awk '{print $1"\t"$3}' S1.covs.tsv | grep -v '^#' > S1.abundance.txt

### second: run maxbin2
conda activate py36
run_MaxBin.pl -thread 4 -contig S1.contigs.fa -abund S1.abundance.txt -out S1/S1

### CheckM

# (out).0XX.fasta	the XX bin. XX are numbers, e.g. out.001.fasta
# (out).summary	summary file describing which contigs are being classified into which bin.
# (out).log	log file recording the core steps of MaxBin algorithm
# (out).marker	marker gene presence numbers for each bin. This table is ready to be plotted by R or other 3rd-party software.
# (out).marker.pdf	visualization of the marker gene presence numbers using R
# (out).noclass	all sequences that pass the minimum length threshold but are not classified successfully.
# (out).tooshort	all sequences that do not meet the minimum length threshold.

##################   dbcan   ################
### 
# time_consuming: 
cd ..
mkdir 07_dbcan
# use protein sequence to find CGCs
ln -s ../03_prokka/*/*.faa ./
run_dbcan.py S1.faa protein --out_dir S1_dbcan_out/ --db_dir /vd02/home2/Xue/db/ --dia_cpu 10 --hmm_cpu 10 --tf_cpu 10
# use DNA sequence to find CGCs
ln -s ../02_megahit/*/*.fa ./
run_dbcan.py S1.contigs.fa meta --out_dir S1_dbcan_out/ --db_dir /vd02/home2/Xue/db/ --dia_cpu 10 --hmm_cpu 10 --tf_cpu 10
# dbcan apply hmm, diamond, and hotpep methods to predict CGCs. the CGCs predicted by only one methods are not included
# 

### singleM - alpha diversity estimation
ln -s ../0_rawdata/* ./
singlem pipe --forward trimmed.S1_r1.fq.gz --reverse trimmed.S1_r2.fq.gz --otu_table s1.tsv --threads 2 --output_extras

