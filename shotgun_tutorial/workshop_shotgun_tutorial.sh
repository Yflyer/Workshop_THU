### no-thinking workflow for processing single metagenomic data ###
### code author: Yufei Zeng
### yfzeng0827@hotmail.com
### github.com/Yflyer/

##### before the work, please copy the data to your project directory
mkdir [your_project_name]
cp -r /vd04/yufei/workshop_shotgun/test/0_rawdata [your_project_name]

###################  qc  #################
conda activate py36
mkdir 01_cleandata
cd 01_cleandata
ln -s ../0_rawdata/* ./
trimmomatic PE -phred33 -threads 4 S1_r1.fq.gz S1_r2.fq.gz \
      trimmed.S1_r1.fq.gz outtrimmed.S1_r1.fq.gz trimmed.S1_r2.fq.gz outtrimmed.S1_r2.fq.gz  \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 \
      MINLEN:50

##################   contig assembly   #######################
### we used megahit to save time and memory
cd ..
mkdir -p 02_megahit
cd 02_megahit
ln -s ../01_cleandata/trimmed* ./
megahit -1 trimmed.S1_r1.fq.gz -2 trimmed.S1_r1.fq.gz --min-count 2 --k-list 29,39,51,67,85,107,133 -m 0.5 -t 4 --min-contig-len 500 --out-prefix S1 -o S1

sed -i "s/>/>S1\_/1" S1/S1.contigs.fa

#################   ORF prediciton  ############
### we used prokka to avoid too much detail
cd ..
mkdir -p 03_prokka
cd 03_prokka
ln -s ../02_megahit/*/*.fa ./
prokka --metagenome --cpus 4 --addgenes --outdir S1 --prefix S1 --mincontiglen 500 S1.contigs.fa

sed -i "s/>/>S1\_/1" S1/S1.ffn
sed -i "s/>/>S1\_/1" S1/S1.faa

##################   mapping  ################
### we used bamm to quickly get mapping and coverage
cd ..
#switch to conda env py27
conda activate py27
mkdir -p 04_mapping
cd 04_mapping
ln -s ../03_prokka/*/*.ffn ./
ln -s ../01_cleandata/trimmed* ./

bamm make -d S1.ffn -c trimmed.S1_r1.fq.gz trimmed.S1_r2.fq.gz -t 4 --out_folder S1
bamm parse -c S1.covs.tsv -l S1.links.tsv -b S1/S1.ffn.trimmed.S1_r1.bam

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

################################################
### we used mmseq to quickly cluster
cd ..
conda activate py36
mkdir 05_orf_clustering
cd 05_orf_clustering

cat ../03_prokka/*/*.faa >> merge.faa

mkdir db clu_db clu_rep
mmseqs createdb merge.faa db/faa
mmseqs linclust db/faa clu_db/faa tmp --threads 4 --min-seq-id 0.8
mmseqs createtsv db/faa db/faa clu_db/faa clu.tsv --threads 4

mmseqs createseqfiledb db/faa clu_db/faa clu_rep/faa
mmseqs result2flat db/faa db/faa clu_rep/faa clu_rep.fasta
#################################################

### and next, we change to R environment to process our ORFs data