# Simplified Tutorial of Processing Sequencing Data

Yufei Zeng
yfzeng0827@hotmail.com
github.com/Yflyer/
THU， 15/03/2020

-----


***In this turtorial, we will operate the sequecning data of microbial community, to generate the relative abundunce table of species and their taxanomy information.***

### 1. Linux system

###### 1.1 Information for your sign-in

Please,  
(1) connect to THU vpn in advance  
(2) install a remote terminal (e.g., Xshell)

```
IP：166.111.42.42
User: test0
password: ********
```
*Why we use Linux?*
* High stability
* Open source
* Perfect For Programing

###### 1.2 Basic Commands of Linux

<br> Commands always look like: <br>
`Command + [space] + [option/ -* / --***] + [space] + [objects: file/ parameters/ link/ path]
`

---

|Command|Description|common option|
|-|-|-|
|cd|Enter the directory| |
|ls|View directory contents |-a, -l|
|mkdir|Make a new directory||
|rm|Delete|-r, -f
|vim |Edit files with vim|
|cat|Merge or view files|-b,
|head|View file header|-n|
|mv|Cut & Paste (Rename)||
|pwd| Show current path||
|which|View the command execution path||
|locate|Indexing Files or Folders||
|**man**|View the command documentation||

---
\# The path in Linux

\# relative path  
From the current directory: `./file_name  `  
Upper directory: `../file_name`

\# Absolute path  
From the root directory /: `/home/test`  
Home directory shortcuts: `~`

---

\# Example: view your fastq data
```
# Enter the directory of forward sequence (R1)
cd R1

# Extract sample list
ls | grep 18S

# Extract 4 lines per sample
head -4 $(ls | grep 18S) > data_report.txt

# View data report
cat data_report.txt

# optional:directly view the data by less  
less 18S_R1_X2016_10N
```
Two common type of sequence data format:

\# **fastq**
```
@M01056:178:000000000-AKTR1:1:1101:17173:1474 1:N:0:0
TACATAGTGCCAGCAG…
+
3>>AAFFFFFFFGGGGG…
```
\# **fasta**
```
>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMA
```


---

###### 1.3 Use conda

Conda is an open-source package management system and environment management system   that runs on Linux.
   * Conda quickly installs, runs, and updates packages and their dependencies

   * Conda easily creates,saves, loads, and switches between environments on your  local computer

   * Conda can package and distribute software for any language and any version of them

\# common conda mangement command
   ```
   # show current information of conda
   conda info

   # list all installed environment
   conda env list
   ```

\# **Optional: Configure QIIME2 by conda**
```
# download the yml
 wget https://data.qiime2.org/distro/core/qiime2-2020.2-py36-linux-conda.yml

# configure env by yml (time-consuming!!!)
 conda env create -n qiime2-2020.2 --file qiime2-2020.2-py36-linux-conda.yml

# activate QIIME2 environments
 conda activate qiime2-2020.2
```
### 2. QIIME2
###### *About QIIME2*
QIIME2 is a next-generation microbiome bioinformatics platform that is extensible, free, open source, and community developed.

Official Tutorial: https://docs.qiime2.org/2019.10/tutorials/
中文版教程： https://blog.csdn.net/woodcorpse/article/details/77929607

###### 2.1 Import Data

   * data type: split files of fastq
   * data source: The University of Oklahoma
   * method: 18s rRNA high through-put amplicon sequencing
   * strategy: pair-end 250 bp without primer on V4 region (**(F565-R981, 416 bp)**)
   * forward primer: CCAGCASCYGCGGTAATTCC **(20 bp)**
   * reverse primer: ACTTTCGTTCTTGATYRA **(18 bp)**
   * The pre-required files: R1 data, R2 data, and mapping information

```
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2
```

\# Preview the sequence quality
```
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
```
A easy way to preview it: https://view.qiime2.org/

###### 2.2 Denoise data by DADA2

The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used.
```
#evalutate the length of sequence in bad quality
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 43 \
  --p-trim-left-r 42 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 200 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```  
* `--p-trim-left-f` or `--p-trim-left-r`: the length for trimming beginning part of forward/reverse sequence  
* `--p-trunc-len-f` or `--p-trunc-len-r`:  the length for cut off end part of forward/reverse sequence  

![avatar](Sequence.png)

*Think a question:*  
*The length of target fragment are 981 - 565 - (20+18) = 378 bp*  
*So what is the length of the merged sequence after denoising in this step?*

\# optional: visualize the DADA2 result
```
# OTU visualization
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

# representative sequences visualization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# denoising information
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```
\# optional: export the biom file
```
# export biom file
qiime tools export  \
  --input-path table.qza \
  --output-path result
# convert biom to tsv
biom convert -i result/feature-table.biom -o feature-table.tsv --to-tsv
# export representative sequences
qiime tools export  \
  --input-path rep-seqs.qza\
  --output-path result
```
###### 2.3 Taxanomy classifier
We will train the Naive Bayes classifier using silva_132 reference sequences (clustered at 99% similarity) and classify the representative sequences from the Moving Pictures dataset.

```
mkdir training_classifiers
cd training_classifiers
```
\# Two elements are required for training the classifier: the reference sequences and the corresponding taxonomic classifications.
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path silva_132_99_18S.fna \
  --output-path 99_silva_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path silva_taxonomy_7_levels.txt \
  --output-path ref_silva_taxonomy.qza
```
 \# taxonomic classification accuracy of 16S rRNA gene sequences improves when a Naive Bayes classifier is trained on only the region of the target sequences that was sequenced (Werner et al., 2012).
```
qiime feature-classifier extract-reads \
  --i-sequences 99_silva_otus.qza \
  --p-f-primer CCAGCASCYGCGGTAATTCC \
  --p-r-primer ACTTTCGTTCTTGATYRA \
  --p-min-length 300 \
  --p-max-length 450 \
  --o-reads ref_silva_seqs.qza
```
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref_silva_seqs.qza \
  --i-reference-taxonomy ref_silva_taxonomy.qza \
  --o-classifier silva_99_classifier.qza
```
\# classify the representative sequences from previous dataset
```
qiime feature-classifier classify-sklearn \
  --i-classifier  silva_99_classifier.qza \
  --i-reads ../rep-seqs.qza \
  --o-classification taxonomy.qza
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

### Summary
In this tutorial, we achieved following goals:
* Use Linux
* configure environment for processing dataset
* change the unstructural sequence data into structual dataframe
* annotate the taxanomy information of species
