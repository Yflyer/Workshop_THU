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

mkdir binning_data
cd binning_data
cat *.fa > merge.contigs.fa
conda activate py36
seqkit seq -m 1500 merge.contigs.fa > merge.trimlen.fa
ln -s ../01_cleandata/trimmed.* .

conda activate py27
bamm make -d merge.trimlen.fa -c trimmed.10-V1_R1.fq.gz trimmed.10-V1_R2.fq.gz 10-V2_R1.fq.gz trimmed.10-V2_R2.fq.gz trimmed.10-V3_R1.fq.gz trimmed.10-V3_R2.fq.gz trimmed.10-V4_R1.fq.gz trimmed.10-V4_R2.fq.gz trimmed.10-V5_R1.fq.gz trimmed.10-V5_R2.fq.gz -t 60 --out_folder merge
