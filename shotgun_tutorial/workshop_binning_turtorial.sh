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

mkdir **
cd **
ln -s ../02_megahit/1-*/*.fa .
ln -s ../01_cleandata/L1-* .
conda activate py36
seqkit seq -m 1500 *.fa > L1.fa

conda activate py27
bamm make -d L1.fa -c L1-*.fq.gz -t 60 --out_folder L1
bamm parse -c L1.cov.tsv -b L1/*.bam -t 60

# file excluding first N lines: tail -n +<N+1> <filename> 
# awk default print: awk "{}1"
cat L10.cov.tsv | tail -n +2 | awk '{$2=""}1' > L10.abund.tsv

conda activate py36
mkdir L10.binning
run_MaxBin.pl -contig L10.fa -abund L10.abund.tsv -out L10.binning/L10 -thread 60

conda activate py27
bamm make -d L10.fa -c 10-V1.r1.fq.gz 10-V1.r2.fq.gz 10-V2.r1.fq.gz 10-V2.r2.fq.gz 10-V3.r1.fq.gz 10-V3.r2.fq.gz 10-V4.r1.fq.gz 10-V4.r2.fq.gz 10-V5.r1.fq.gz 10-V5.r2.fq.gz -t 60 --out_folder L10







