
# Shotgun data processing
Yufei Zeng
yfzeng0827@hotmail.com
github.com/Yflyer/  
THU， 12/03/2019

-----
## 1. Config the conda
(i) （国内）调整linux apt默认安装源（清华源教程）
Ubuntu 的软件源配置文件是 /etc/apt/sources.list。将系统自带的该文件做个备份，将该文件替换为下面内容，即可使用 TUNA 的软件源镜像。
```
# 默认注释了源码镜像以提高 apt update 速度，如有需要可自行取消注释
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic-updates main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic-updates main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic-backports main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic-backports main restricted universe multiverse
deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic-security main restricted universe multiverse
# deb-src https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ bionic-security main restricted universe multiverse
```
(ii) （国内）调整conda默认安装源（清华源教程）
创建修改用户目录下的 .condarc 文件
```
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  msys2: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
```
optional: 设置包安装时候的源可见
```
conda config --set show_channel_urls yes
```
重启terminal，检查conda配置是否完成
```
conda config --show
```


## 1. Environment Preparation
创建一个python3.6环境避免过多冲突
```
conda create -n py36 python=3.5
```

```
# QC software  
conda install -n py36 -c bioconda fastqc
conda install -n py36 -c bioconda trimmomatic
conda install -n py36 -c bioconda bbmap
conda install -n py36 -c bioconda multiqc
conda install -n py36 -c bioconda khmer
# assembly software  
conda install -n py36 -c bioconda megahit
conda install -n py36 -c bioconda spades
```
* 看到一个很奇葩的[bug](http://www.mamicode.com/info-detail-2272598.html)，如果没有桌面是无法运行fastqc交互界面，server版的linux还要安装桌面（图形界面）。ORZ...我忙了好久，发现这个bug根本不影响，只是打不开fastqc图形版而已......
```
# update sources list
sudo apt-get update
# install
sudo apt-get install gnome
# optional if interuption
sudo apt-get install --reinstall gnome
```
安装JGI的bbtools，如果版本冲突，可以用anaconda search寻找合适版本
```
# search compatiable bbtools  

# install bbtools from JGI
conda install -c agbiome bbtools
```
## 2. Quality Control
fastqc  
trimmomatic  
adapter sequence files包括seq2和seq3，使用最新的seq3版本  
multiqc  

## 3. *De novo* Assembly

### *About Linux*
find
grep

## 4. toolkit env config
### kingfisher
conda create -c conda-forge -c bioconda -n kingfisher pigz python extern curl sra-tools pandas requests aria2conda activate kingfisher

\#使用conda activate不能成功激活环境时可以尝试使用：
\# source activate kingfisherpip install bird_tool_utils'>='0.2.17git clone https://github.com/wwood/kingfisher-downloadcd kingfisher-download/binexport PATH=$PWD:$PATHkingfisher -h#弹出帮助文档即安装成功

### das_tool
conda create -c conda-forge -c bioconda -n das_tool das_tool

### dbcan
#### dbcan database config
cd db \
    && wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312019.fa.nr && diamond makedb --in CAZyDB.07312019.fa.nr -d CAZy \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt && mv dbCAN-HMMdb-V8.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm \
    && cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
    && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff

# Check Program.
run_dbcan.py EscheriaColiK12MG1655.fna prok --out_dir test –db_dir /vd02/home2/Xue/db/

### MetaBAT2
conda activate py27

###  checkM-genome
conda activate py36
pip3 install numpy
pip3 install matplotlib
pip3 install pysam
pip3 install checkm-genome
