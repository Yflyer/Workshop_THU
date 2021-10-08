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

# Banfield lab:
# BBmap remove Illumina adaptor and phiX sequences
# Sickle to qc control
# individually assembled by IDBA-UD
mkdir L1
cd L1
ln -s ../02_megahit/1-*/*.fa .
ln -s ../01_cleandata/L1-* .

### individual assembly stratgy
conda activate py36
seqkit seq -m 1500 *.fa > L1.fa

conda activate py27
bamm make -d L1.fa -c L1-*.fq.gz -t 60 --out_folder L1
bamm parse -c L1.cov.tsv -b L1/*.bam -t 6

# file excluding first N lines: tail -n +<N+1> <filename>
# awk default print: awk "{}1"
cat L1.cov.tsv | tail -n +2 | awk '{$2=""}1' > L1.abund.tsv

conda activate py36
mkdir L1.binning
run_MaxBin.pl -contig L1.fa -abund L1.abund.tsv -out L1.binning/L1 -thread 6

### check M
checkm lineage_wf <bin folder> <output folder>

### co-assembly startgy by SPAdes
mkdir L1.spades
cd L1.spades
# merge at first (the same as gzip)
pigz -d -c  ../*1.fq.gz | pigz -c > merged.R1.fq.gz
pigz -d -c  ../*2.fq.gz | pigz -c > merged.R2.fq.gz

### take a try
### bbnorm
### bbnorm to adjust coverage will be good for spades:
### http://seqanswers.com/forums/showthread.php?t=49763
# loglog.sh in1=merged.R1.fq in2=merged.R2.fq # ~60G 5000s
bbnorm.sh in1=merged.R1.fq in2=merged.R2.fq out=highpass.fq outt=lowpass.fq passes=1 target=40 min=5 -Xmx250g threads=40 # consume 230G memory actually, ~3 hours
# after bbnorm, files been interleaved
bbmerge-auto.sh in=highpass.fq out=bbmerge.fq outu=unmerged.fq rem extend2=50 k=62 -Xmx250g -t=40
spades.py --meta --12 bbmerge.fq -o L1 --threads 8 -m 250

###
mkdir L1
ln -s ../L1-* .
ls -d *.R1.fq.gz | cut -d '.' -f1 > L1_list.txt

### bbduk
# used to remove contaminants, reads that contained adapter sequences, and right quality trim reads where quality dropped to 0.
# BBDuk was also applied to eliminate reads containing 1 or more “N” bases, having an average quality score across the read of less than 13 or containing a minimum length of ≤ 41 bp or 33% of the full read length.
# Using BBMap, reads that were mapped to masked human, cat, dog, and mouse references at 95% identity and aligned to common microbial contaminants were separated.
# right end, Kmer length = 27, maximum substitutions = 1, minimum quality = 20, minimum overlap = 20, minimum length = 20
### bbduk + bbmap
parallel -j 11 --xapply 'bbduk.sh in1={1}.R1.fq.gz in2={1}.R2.fq.gz out={1}.rmadp.fq ref=~/bbmap_resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 -t=6 tpe tbo' :::: L1_list.txt
parallel -j 11 --xapply 'bbduk.sh in={1}.rmadp.fq out={1}.qc.fq qtrim=r trimq=20 maq=20 -t=6' :::: L1_list.txt
parallel -j 11 --xapply 'bbduk.sh in={1}.qc.fq out={1}.clean.fq ref=~/bbmap_resources/phix.fa k=31 hdist=1 stats={1}.stats.txt -t=6' :::: L1_list.txt
#To index and map at the same time:
# To build an index in-memory without writing to disk:
# bbmap.sh in=${i}.qc.fq out=${i}.mapped.sam ref=ref.fa nodisk -t=40
# To split input into mapped and unmapped, in fastq format:
# bbmap.sh in=${i}.mapped.sam outu=u${i}.rmhos.fq -t=40

bbduk.sh in=reads.fq out=unmatched.fq outm=matched.fq ref=phix.fa k=31 hdist=1 stats=stats.txt

### bbmerge
# It is designed for kmer-based operations using Tadpole, which include both merging overlapping and non-overlapping reads, kmer-based error-correction, and kmer-based filtering.
# Kmer-based operations should only be used with shotgun (randomly-fragmented) libraries, never with amplicon libraries (such as 16S).
#  If you run BBMerge, and under, say, 15% of the reads merge, even at very loose stringency, it’s probably a waste of time to merge – you’ll just make the workflow more complicated, and possibly get a lot of false-positives. Also, don’t try to merge single-ended libraries or long-mate-pair libraries that are not in an “innie” orientation.


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
