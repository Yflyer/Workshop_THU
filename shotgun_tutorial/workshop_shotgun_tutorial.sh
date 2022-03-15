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

#### before your work #######
# prepare necessary materials:
# Copy bbmap_resources to your home path

#### activate conda env #####
conda activate py36
###################  qc  ######################
# time_consuming: ★☆
mkdir 01_cleandata
cd 01_cleandata
ln -s ../0_rawdata/* ./ # link
bbduk.sh in1=S1_r1.fq.gz in2=S1_r2.fq.gz out=S1.clean.fq ref=~/bbmap_resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=r trimq=20 maq=20 -t=4 tpe tbo

##################   contig assembly   #######################
### we used megahit to save time and memory
# time_consuming: ★★★★
cd ..
mkdir -p 02_megahit
cd 02_megahit
ln -s ../01_cleandata/*clean.fq ./
megahit --12 S1.clean.fq --min-count 2 --k-list 29,39,51,67,85,107,133 -m 0.5 -t 4 --min-contig-len 500 --out-prefix S1 -o S1

# add prefix to each contig by sed 
# this step is very important to make contig ID unique across samples
sed -i "s/>/>S1\_/1" S1/S1.contigs.fa

#################   ORF prediciton by prokka #######################
### we used prokka to avoid too much detail
# time_consuming: ★★☆
cd ..
mkdir -p 03_prokka
cd 03_prokka
ln -s ../02_megahit/*/*.fa ./
prokka --metagenome --cpus 4 --outdir S1 --prefix S1 --mincontiglen 500 --locustag S1 S1.contigs.fa


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

##################  ORF: Cluster and Mapping  ################
### we used mmseq to quickly cluster
# time_consuming: ★★☆
cd ..

mkdir -p 04_ORF_mapping
cd 04_ORF_mapping
ln -s ../01_cleandata/*clean.fq ./
ln -s ../03_prokka/*/*.ffn ./
ln -s ../03_prokka/*/*.ffa ./

### cluster ORF
mkdir db clu_db clu_rep
# faa cluster (to get cluster seqs of ffn in mapping)
# We set the similiarity to 0.9
mmseqs createdb *.faa db/faa
mmseqs linclust --min-seq-id 0.9 db/faa clu_db/faa tmp --threads 24
mmseqs createtsv db/faa db/faa clu_db/faa clu.faa.tsv --threads 24 # we could use clu.faa.tsv to extract all seq we need by seqkit

### get representative sequences
mmseqs createsubdb clu_db/faa db/faa clu_rep/faa
mmseqs convert2fasta clu_rep/faa clu_rep.faa

### extract rep fnn by faa ID
seqkit seq -j 24 clu_rep.faa -n -i > rep_list.txt
seqkit grep -j 24 --pattern-file rep_list.txt *.ffn > clu_rep.ffn

### we used coverm to quickly get mapping and coverage
# switch to conda env coverM: calculate coverage table
# time_consuming: ★★★★
conda activate coverM
coverm contig --interleaved *.fq --reference clu_rep.ffn -t 60 --bam-file-cache-directory orf_bam -o clu_rep_coverage.tsv

### The table you get might need extra qulaity control includ but not limit to singleton removal. Please do it in you R

##################  ORF Annotation   ################
### There are a lot of annotation methods. We take cazymes and KEGG as example.
### !NOTE: IF your clu_rep.faa is too LARGE, please use seqkit to split them into small files and then annotation
### !NOTE: You can use parallel to annotatate small files quickly!
### we used mmseq to quickly cluster
# time_consuming: ★★★
cd ..
mkdir 05_ORF_annotation
cd 05_ORF_annotation
### cazyme annotation
ln -s ../04_ORF_mapping/clu_rep.* ./

conda activate run_dbcan

run_dbcan.py clu_rep.faa protein --out_pre clu_rep.faa --out_dir cazyme_annoation --db_dir /vd03/home/MetaDatabase/dbcan --dia_cpu 1 --hmm_cpu 1 --tf_cpu 1 --hotpep_cpu 1

# check the result: overview.txt

### KEGG annotation 
# time_consuming: ★★★★
conda activate py36
exec_annotation -f  detail-tsv -E 1e-5 --profile /vd03/home/MetaDatabase/KOfam_2019/Kofam/profiles/ --ko-list /vd03/home/MetaDatabase/KOfam_2019/Kofam/ko_list --cpu 48 --tmp-dir ./ko_tmp -o KEGG_annoation/clu_rep.txt clu_rep.ffn
# check the result: clu_rep.txt

