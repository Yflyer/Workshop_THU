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
# maxbin
mkdir L1.maxbin
run_MaxBin.pl -contig L1.fa -abund L1.abund.tsv -out L1.maxbin/L1 -thread 6

# metabat
conda activate metabat2

runMetaBat.sh -t 20 -m 1500 -o L1.metabat/L1 L1/L1.fa L1/L1/*.bam
metabat2 -i L1/L1.fa -a L1.fa.depth.txt -t 6 -m 1500 -o L1.metabat/L1
runMetaBat.sh -t 20 -m 1500 -o ${i}.metabat/${i} ${i}/${i}.fa

ls -ad L* > site_list.txt
while read i; do
mkdir ${i}.metabat
jgi_summarize_bam_contig_depths --outputDepth ${i}.fa.depth.txt ${i}/${i}/*.bam
metabat2 -i ${i}/${i}.fa -a ${i}.fa.depth.txt -t 20 -m 1500 -o ${i}.metabat/${i}
done <site_list.txt

#
# concot
cut_up_fasta.py L1/L1.fa -c 10000 -o 0 --merge_last -b L1_10K.bed > L1_10K.fa
concoct_coverage_table.py L1_10K.bed L1/L1/*.bam > L1_concoct_table.tsv
concoct --composition_file L1_10K.fa --coverage_file L1_concoct_table.tsv -b L1_concoct/
merge_cutup_clustering.py L1_concoct/L1_clustering_gt1000.csv > L1_concoct/L1_clustering_merged.csv
mkdir L1_concoct/fasta_bins
extract_fasta_bins.py L1/L1.fa L1_concoct/L1_clustering_merged.csv --output_path L1_concoct/fasta_bins
# ESOM: waiting for update, because three methods is sufficient for das-tool comparison

# metawrap
metawrap binning -o ${i}.metawrap -t 20 -a ${i}/${i}.contigs.fa --metabat2 --maxbin2 --concoct --interleaved ${i}/*.fastq


### das_genome
conda activate das_tool
Fasta_to_Scaffolds2Bin.sh -i L1.maxbin/ -e fasta > L1.maxbin2das.tsv
Fasta_to_Scaffolds2Bin.sh -i L1.metabat/ -e fasta > L1.metabat2das.tsv
parallel -j 4 --xapply 'Fasta_to_Scaffolds2Bin.sh -i {1}.metabat -e fa > {1}.metabat2das.tsv' :::: site_list.txt
parallel -j 4 --xapply 'Fasta_to_Scaffolds2Bin.sh -i {1}.maxbin -e fasta > {1}.maxbin2das.tsv' :::: site_list.txt

DAS_Tool -i L1.metabat2das.tsv,L1.maxbin2d as.tsv -l metabat,maxbin -c L1/L1.fa -o L1.dastool -t 15

# LOOP
ls -d L* | cut -d '.' -f1 > site_list.txt
while read i; do
  mkdir ${i}.dastool
  DAS_Tool -i ${i}.metabat2das.tsv,${i}.maxbin2das.tsv -l metabat,maxbin -c ${i}/${i}.fa -o ${i}.dastool/${i} -t 12 --write_bins 1
done <site_list.txt

### check M
# pplacer_threads: number of threads used by pplacer (memory usage increases)
# nt: generate nucleotide gene sequences for each bin
checkm lineage_wf --nt -x fa -t 8 --pplacer_threads 4 -f test.tsv --tab_table L1.dastool/L1_DASTool_bins L1.checkm
### check M for CPR:
# Identify marker genes in bins and calculate genome statistics
checkm analyze -x fa -t 8 ~/hmm_sets/cpr_43_markers.hmm L1.dastool/L1_DASTool_bins L1.checkm.cpr
# Assess bins for contamination and completeness.
checkm qa -t 8 -o 1 -f test.cpr.tsv --tab_table ~/hmm_sets/cpr_43_markers.hmm L1.checkm.cpr

# LOOP
while read i; do
  checkm lineage_wf --nt -x fa -t 8 --pplacer_threads 4 -f ${i}.check.tsv --tab_table ${i}.dastool/{i}_DASTool_bins ${i}.checkm
done <site_list.txt

### mimag selection
# minimum information about a metagenomeassembled genome (MIMAG) standards:
# high: >90% completeness and <5% contamination, presence of 5S, 16S and 23S rRNA genes, and at least 18 tRNAs;
# medium: ≥ 50% completeness and <10% contamination.
# Ref: Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life
# rename file and seqs in dastool output
while read i; do
  #echo *.dastool/*_DASTool_bins/${i}.fa
  #sed -i "s/>/>${i}\_/1" ${i}.dastool/${i}_DASTool_bins/*.fa
  cd ${i}.dastool/${i}_DASTool_bins
  for m in *.fa
  do
    mv ${m} ${i}_${m}
  done
  cd ../..
  #mv ${i}.dastool/${i}_DASTool_bins/*.fa ${i}.dastool/${i}_DASTool_bins/${i}_*.fa
done <site_list.txt

# process checkM result by awk command
touch MIMAG.tsv
while read i; do
  awk -v var=${i} '$13 > 50 && $14 <10 {printf "%s\t%s\t%s\t%s\t%s\n",var,var"_"$1,$2,$13,$14}' ${i}.check.tsv >> MIMAG.tsv
done <site_list.txt

# get MIMAG by checkM from DAStool
cat MIMAG.tsv | cut -f 1-2 > MIMAG_list.txt
awk '{print $1".dastool/"$1"_DASTool_bins/"$2".fa";}' MIMAG.tsv > MIMAG_list.txt
awk '{print $2;}' MIMAG.tsv | cut -d '_' -f1,2 --output-delimiter=' ' | awk '{print $1".checkm/bins/"$2;}' > MIMAG_checkm.list

### Taxanomy annotation of binning genome: tools, biomes, data features, reference, and notes.
# GTDB  soils  bacterial binned metagenome  EM: similar functional profiles of CPR phyla in soils  2020 citation 677
# Phylosift  hydrothermal_sediments prokayotic binned metagenome  NC: expansive metabolic versatility in GB hydrothermal sediments  Only phylogeny?
# mOTUs marker  Ocean  read-based genes mapping  Science: Structure and function of the globalocean microbiome 2019 citation 123
# Metaxa2  soils read- or contig- based 16s NM: Candidatus Udaeobacter copiosus

### taxa anotation: GTDB ###
mkdir 07_taxa 
cd 07_taxa
ln -s ../06_bin/MIMAG_list.txt .
# cp MIMAG
mkdir gtdb
while read i; do
  cp ../06_bin/${i} gtdb
done <MIMAG_list.txt
conda activate gtdbtk-1.5.0
gtdbtk classify_wf --cpus 60 -x fa --pplacer_cpus 1 --genome_dir gtdb --out_dir gtdb_result
# get mimag stat by bbmap
stats.sh in=$(ls *.fa | paste -s -d ",") format=6 threads=16 addname=t
statswrapper.sh in=$(ls *.fa | paste -s -d ",") format=6 threads=16 addname=t gc=mimag.gc.txt gcformat=4

# build tree
# followed by the ‘phylosift align’ mode. The concatenated protein alignments of 37 elite marker genes (concat.updated.1.fasta) were combined for all genomes of interest and trimmed using TrimAL (version 1.2) using the automated1 setting65. A phylogenetic tree was generated using a maximum likelihood-based approach using RAxML (version 8.2.10, called as: raxmlHPC-PTHREADS-AVX -f a -m PROTGAMMAAUTO -N autoMRE)
# The resulting alignments were stripped of columns containing >95% gap positions. Individual stripped alignments were concatenated and a phylogenetic tree was constructed using RAxML v8.2.1076 on the CIPRES Science Gateway{Miller:vv}. RAxML was called as follows: raxmlHPC-HYBRID -s input -N autoMRE -n result -f a -p 12345 -x 12345 -m PROTCATLG.
conda activate py36
raxmlHPC -s gtdbtk.bac120.user_msa.fasta -T 60 -N autoMRE -n test -f a -p 12345 -x 12345 -m PROTCATLG

cd ..
### cazy annotation

mkdir 08_dbcan
cd 08_dbcan
# use protein sequence to find CGCs
ln -s ../07_taxa/mimag .
ln -s ../06_bin/MIMAG_list.txt .

cd mimag
conda activate run_dbcan
for i in *.fa
do
  run_dbcan.py ${i} meta --out_dir ${i/\.fa/} --db_dir /vd03/home/MetaDatabase/dbcan --dia_cpu 16 --hmm_cpu 16 --tf_cpu 16 --hotpep_cpu 16
done

ls -d * > fna_list.txt
awk 'NR>1 {print $0;}' overview.txt | awk '$5 > 2 {print $2;}' | cut -d '(' -f1 | cut -d '_' -f1 | paste -s -d ';'

touch cazy_result.tsv
while read i; do
  gene=$(awk 'NR>1 {print $0;}' ${i}/overview.txt | awk '$5 > 2 {print $2;}' | cut -d '(' -f1 | cut -d '_' -f1 | paste -s -d ';')
  echo $i $gene >> cazy_result.tsv
done <fna_list.txt


#################  KOfam  ################
#在国家微生物科学数据中心网站可以下载最新版2019年KOfam ko_list和profiles （KEGG Orthologs（KOs）的定制HMM数据库）为用户的序列数据的搜索，通过将用户的序列数据与KEGG路径和EC编号联系起来，得到注释结果。
mkdir orf
conda activate py36
ln -s ../../06_bin/MIMAG_checkm.list
while read i; do
  cp -r ../../06_bin/${i}* .
done <MIMAG_checkm.list

### 48 core, 24 hour, complete nearly 80 binned genome; slowly
for i in $(ls -d *); do
    exec_annotation -f  detail-tsv -E 1e-5 --profile /vd03/home/MetaDatabase/KOfam_2019/Kofam/profiles/ --ko-list /vd03/home/MetaDatabase/KOfam_2019/Kofam/ko_list --cpu 48 --tmp-dir ./ko_tmp -o ${i}_kofam.txt ${i}/genes.fna
done

-v var=${i} $13 > 50 && $14 <10 
awk '{printf $3}' L2.52_kofam.txt
sed -i "s/K//1" L2.52_kofam.txt
sed -n '/K00174/p' L2.52_kofam.txt

K00174
K00169
K00123
K00128
K00001
K00925
K01905

K00399
K00198
K14138
K15023
K00192
K00195


### genome coverage
cat mimag/*.fa > totalMAG.fa
mkdir reads
ln -s ../01_bbmap/*fq reads/
conda activate coverM

# Dereplicate at 99% (after pre-clustering at 95%) a directory of .fna
coverm cluster -x fa -t 60 --genome-fasta-directory mimag --output-representative-fasta-directory drep_mimag --output-cluster-definition mimag_clusters.tsv
# cal cov
coverm genome -t 60 --interleaved reads/*fq -x fa --genome-fasta-directory mimag --bam-file-cache-directory mimag_bam -o MAG_coverage.tsv

### iTOL
# label: tree name

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