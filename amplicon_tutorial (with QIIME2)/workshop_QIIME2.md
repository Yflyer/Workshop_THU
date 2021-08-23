
# Linux, QIIME2, & DADA2
Yufei Zeng
yfzeng0827@hotmail.com
github.com/Yflyer/
THU， 11/15/2019

-----
12/02/2019 添加了classifier训练的方法和注释教程
12/11/2019 添加了phylogenetic inference的方法教程
12/16/2019 添加了Alpha/Beta diversity的计算代码

## 1. Install your linux
For win10 user: enable the WSL and install Ubuntu 18.04 from Microsoft Store
https://docs.microsoft.com/en-au/windows/wsl/install-manual
（虽然说了有点多余）装系统后新建用户输入密码时，是不会显示密码的，不要以为出了bug，不要关闭或者跳过，输错了再输就是了。
For MacOS user: almost nothing to do
### *About Linux*
1991年，**由于贫穷**，初生牛犊的北欧小伙Linus单枪匹马一个人基于Unix系统编写了Linux系统 —— 一个完全免费开源的高效开发系统，对整个计算机行业产生了巨大而深远的影响。
* Unix系统作为多任务系统的祖师爷级产品，由贝尔实验室的geeker们研发，逐渐成为行业标准，后来苹果的MacOS系统也基于该系统研发（MacOS几乎能和Linux无缝衔接的重要原因）
* 最早的服务器（计算机），比如IBM的701，是无法进行多任务并行操作的，科学家们需要排队使用。随着万维网和Unix的出现，人们可以同时调用多个终端（terminal）去使用服务器。  



## 2. Install Anaconda or Miniconda  

\# optional: show current path
<br>`pwd`<br>

\# download conda<br>
`wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh`
<br> or <br>
`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`  
conda有python2和python3的版本，推荐python3版本的  
*wget命令添加 --O 选项，可以更改文件存放路径及重命名，如wget --O /mnt/e/Anaconda3。sh

\# install conda  
`bash Anaconda3-2019.10-Linux-x86_64.sh`  
or  
`bash Miniconda3-latest-Linux-x86_64.sh`  
* 如果conda.sh的文件是直接拷贝的，上述命令文件名加上路径信息即可，如/mnt/e/Anaconda3-2019.10-Linux-x86_64.sh  
* wsl用户可能需要添加sudo前置，sudo即为管理员权限

\# optional: chown the conda file  
wsl用户sudo后可能需要conda文件夹权限，根据debug提示输入:  
`sudo chown $UID:$GID ~/.conda`  

\# restart the terminal to use conda

\# common conda mangement command
```
conda info #显示conda版本，env path有关的信息
conda update conda #conda 升级  
rm -rf anaconda  # conda卸载（在安装目录下）  
conda env list # 列出所有conda环境  
```  

## 3. Configure QIIME2
\# download the yml  
`wget https://data.qiime2.org/distro/core/qiime2-2019.10-py36-linux-conda.yml `  

\# configure env by yml (time-consuming!!!)  
`conda env create -n qiime2-2019.10 --file qiime2-2019.10-py36-linux-conda.yml
`

\# activate QIIME2 env  
`conda activate qiime2-2019.10`  
*qiime2装的太慢了，我们来了解一下linux吧

### *Let's practice Linux*
在前面短暂的使用中，大家已经接触到了pwd，cd，wget以及安装conda后的conda有关命令，Linux还有各种各样的命令，帮助你即时在命令行界面下，也能像使用windows的图形界面一样完成各种操作。接下来，大家可以尝试在练习这些常用命令的时候，自动匹配在windows界面下的鼠标操作，方便大家记忆。  

\# common command  
Command + [space] + [option/ -* / --***] + [space] + [objects: file/ parameters/ link/ path]  

---

|Command|Description|common option|
|-|-|-|
|cd|进入目录| |
|ls|查看目录内容 |-a, -l|
|mkdir|新建文件夹||
|rm|删除|-r, -f
|vim |用vim编辑文件|
|cat|合并或查看文件|-b,   
|head|查看文件头部|-n|
|mv|剪切黏贴（重命名）||
|man|查看命令说明文档||
|which|查看命令执行路径||
|locate|索引文件或文件夹||

---
\# understand the path in Linux  
\# 相对路径  
从当前目录开始: ./file_name    
上一级目录: ../file_name2：指上一级目录下的文件  
\# 绝对路径  
从根目录/开始: /home/test  
Home目录快捷键：~  

---
## 4. Use QIIME2  
Official Tutorial: https://docs.qiime2.org/2019.10/tutorials/  
中文版教程： https://blog.csdn.net/woodcorpse/article/details/77929607

### 4.1 Data importing
QIIME2支持多种类型的数据导入，比如说单端或者双端、fasta或fastq。由于：i) 目前主流都是双端测序; ii）大部分情况下拿到手都是原始fastq; iii）fasta数据不含质控，无法进行QIIME2中DADA2的聚类，本教程只解决大家导入双端fastq数据的问题。  

**a.后缀名为.fastq.gz的单个压缩文件**  
EMP标准格式，需要有barcode文件  

\# 用gnuzip head等命令查看EMP数据格式的特点  
```
# 进入序列文件夹
cd ~/Workshop_THU/Tutorials_rawdata/P3/emp-paired-end-sequences  
# 解压fastq.gz  
gunzip forward.fastq.gz reverse.fastq.gz barcodes.fastq.gz  
# 每个文件提取4行
head -4 barcodes.fastq forward.fastq reverse.fastq > data_report.txt  
# 查看数据特征
cat data_report.txt
```
* 上述代码可以用vim写入后bash一次执行~   

上述结果可以看到序列之间是没有样品信息的，因此这种格式数据后续处理的时候，需要一个文件提供每个序列所属样本的索引信息，这种文件在QIIME2流程中称为metadata文件，后缀为.tsv

\# use cat, more, or less to preview the metadata!

EMP标准格式的导入代码：  
single-end：
```
qiime tools import \
  --type EMPSingleEndSequences \
  --input-path emp-single-end-sequences \
  --output-path emp-single-end-sequences.qza
```  
paired-end:  
```
qiime tools import \
  --type EMPPairedEndSequences \
  --input-path emp-paired-end-sequences \
  --output-path emp-paired-end-sequences.qza
```  

**b. fastq or fastq.gz的多个文件(Fastq manifest)**  
无barcode文件，文件已经按照样本分成单个小文件(单端或双端)
```
cd ~/Workshop_THU/test_data/R1
ls -a
head -4 *
```
写了一个根据split样品文件夹自动生成R1，R2的manifest的python脚本。在工作文件下输入(3个参数，正向文件夹R1，反向文件R2，输出文件的名字)：  
`python create_manifest.py R1 R2 manifest.tsv`  
可以预览一下效果。  
`cat manifest.tsv`  

paired-end导入代码：
```
# type即Fastq manifest数据
# format一共有四种，根据阈值分类。OK下机的都是33V2
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.tsv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2
```

**c. fastq的单个文件，带tag（如周老师warming数据）**  
该数据fastq首行尾部打入了“--sample_name”，无barcode

\# 对于这一类型的数据，可以先提前在pipeline中切除非生物序列（space_sequence + primer_sequence）,再根据“--sample_name” spilt成单独的文件后导入qiime2。如果觉得使用galaxy pipeline上下传文件麻烦的话，这里也提供了一个split的python脚本（两个参数，第一个为输入文件名，第二个输出文件夹名）,启动后会让你输入用来判断sample name的标识符，此处为“--”（见上）。举个例子：
```
python tag_split.py ***.fastq output_dir
```

* for more information, https://docs.qiime2.org/2019.10/tutorials/importing/

## 4.2 DADA2  
QIIME2中dada2分为paired_end和single_end，此处只展示paired_end处理过程

**a. preview the sequences quality**
```
# 生成qiime2 可视化文件
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
```
* 根据Fastqc的结果在下一步中切割序列
* visualization: https://view.qiime2.org/

**b. set the parameters of DADA2**  
使用的的数据18sv4 (F547-R952),正反序列
```
#evalutate the length of
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 28 \
  --p-trim-left-r 45 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```
\# optional: visualize the DADA2 result
```  
# OTU visualization
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \

# representative sequences visualization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

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


### More information
plugin: https://docs.qiime2.org/2019.10/plugins/

## Summary of using DADA2
* 经过比较，根据序列质量切割后序列保留率明显要高于Btrim后的保留率（~80% vs ~50%）
* pool or not对ASV的结果影响不重要
* DADA2的结果是可合并的

**c. Taxanomy classifier**  
这里选择的ref database是silva的99阈值18s。
```
mkdir training_classifiers
cd training_classifiers
```

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
根据werner et al., 2012 isme.j 上的结果，把ref.seq先裁剪到目标序列后再训练classifier的效果更好。训练方法为传统的naive bayes（galaxy中的RDP的也是这种）。这里我设置18sV4的前后引物，不设置--p-trunc-len。
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
返回上级目录用这个classifier注释样品的ref.seqs（即ASV序列）

```
qiime feature-classifier classify-sklearn \
  --i-classifier  silva_99_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```
## 4.3 Phylogenetic (Diversity) Analysis

建立工作目录
```
mkdir phylogeny
cd phylogeny
```

```
qiime alignment mafft \
  --i-sequences rep-seqs-se.qza \
  --o-alignment aligned-rep-seqs.qza
```  

**b. Alignment**

*(i).masking alignments*
```
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza
```
*(ii).reference alignments*  

PyNAST  
Infernal  
SINA  

**c. Built a tree
FASTTREE  
```
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree fasttree-tree.qza
```
RAxML
```
qiime phylogeny raxml \
  --i-alignment masked-aligned-rep-seqs.qza \
  --p-substitution-model GTRCAT \
  --p-seed 1723 \
  --p-n-searches 5 \
  --o-tree raxml-cat-searches-tree.qza \
  --verbose
```
```
qiime phylogeny iqtree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --p-substitution-model 'GTR+I+G' \
  --o-tree iqt-gtrig-tree.qza \
  --verbose
```

**d. Root the tree
```
qiime phylogeny midpoint-root \
  --i-tree fasttree-tree.qza \
  --o-rooted-tree rooted-fasttree-tree.qza
```
**e. Export the tree
```
qiime tools export \
  --input-path rooted-fasttree-tree.qza \
  --output-path fasttree
```

4.4 Alpha/Beta

qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-fasttree-tree.qza \
    --i-table table.qza \
    --p-sampling-depth 6000 \
    --m-metadata-file manifest.tsv \
    --output-dir core-metrics-results

# beta explore
qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file Manifest.tsv \
    --m-metadata-column Site \
    --o-visualization core-metrics-results/bray_curtis_distance_matrix-site.qzv \
    --p-pairwise

qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file Manifest.tsv \
    --m-metadata-column Status \
    --o-visualization core-metrics-results/bray_curtis_distance_matrix-status.qzv \
    --p-pairwise



### 写在最后
给大家安利一些软件：
文本编辑器： atom，学习成本较高，但是兼容所有语言，拓展性极强
服务器远程软件：Xshell，XFtp，个人版免费使用。
书籍：
鸟哥的LINUX私房菜
