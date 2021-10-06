### no-thinking workflow for processing bining data ###
### code author: Yufei Zeng
### yfzeng0827@hotmail.com
### github.com/Yflyer/

############### before the work ##################
# please copy the data to your project directory
mkdir [your_project_name]
cd [your_project_name]
#cp -r /vd04/yufei/workshop_shotgun/test/0_rawdata [your_project_name]
# make a work screen
# screen is useful, please know more
screen -s binning

# note
# check bamm map detail: orientation insertsize
# check bbmap
# check maxbin2 workflow and parameters
# check metabat2 workflow and example
# check das_tool usage

mkdir L1
cd L1
ln -s ../02_megahit/1-*/*.fa .
ln -s ../01_cleandata/L1-* .

### individual assembly stratgy
conda activate py36
seqkit seq -m 1500 *.fa > L1.fa

conda activate py27
bamm make -d L1.fa -c L1-*.fq.gz -t 60 --out_folder L1
bamm parse -c L1.cov.tsv -b L1/*.bam -t 60

# file excluding first N lines: tail -n +<N+1> <filename> 
# awk default print: awk "{}1"
cat L1.cov.tsv | tail -n +2 | awk '{$2=""}1' > L1.abund.tsv

conda activate py36
mkdir L1.binning
run_MaxBin.pl -contig L1.fa -abund L1.abund.tsv -out L1.binning/L1 -thread 6

### co-assembly startgy by SPAdes
mkdir L1.spades
cd L1.spades
# merge at first (the same as gzip)
pigz -d -c  ../*1.fq.gz | pigz -c > merged.R1.fq.gz
pigz -d -c  ../*2.fq.gz | pigz -c > merged.R2.fq.gz

# try
pigz -d -c -k merged.R1.fq.gz > merged.R1.fq
pigz -d -c -k merged.R2.fq.gz > merged.R2.fq
bbmerge-auto.sh in=reads.fq out=merged.fq outu=unmerged.fq rem extend2=50 k=62
# spades (use less cores for safe!)
spades.py --meta -1 merged.R1.fq.gz -2 merged.R2.fq.gz -o L1 --threads 4 -m 150


# minimum information about a metagenomeassembled genome (MIMAG) standards:
# high: >90% completeness and <5% contamination, presence of 5S, 16S and 23S rRNA genes, and at least 18 tRNAs;
# medium: ≥ 50% completeness and <10% contamination. 
# Ref: Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life


### following procedures:
# rRNA detection: Rfam (INFERNAL V1) cmsearch 1.1.247 (options -Z 1000 --hmmonly --cut_ga)
# tRNA check: tRNAscan-s.e. v.2.049 using the bacterial tRNA model (option -B)

### Genome dereplication:
# all-against-all comparison: MinHash
# clustering: the Mash distance relationships and individual clusters were defined at a cut-off of 0.2
# dereplicating: dRep V2.2

### functional centric target
# For example: WLP in actinobacteria:
# Genome Taxonomy Database (GTDB)-Tk: de_novo_wf --outgroup_taxon p__Chloroflexota --taxa_filter p__Actinobacteria --bac120_ms



