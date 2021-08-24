########
# This R script is written by Yufei Zeng
# yfzeng0827@hotmail.com
# only used for eco925ers
#######
# store the current directory
initial.dir<-getwd()
#########################
input = 'cov'
cluster = 'clu.tsv'


########################
message('=========== =============== (from eco925) ===================== ==============')
message('=========== ========== design for bamm output ================= ==============')
if (!require('dplyr')) {
  message('==== installing required packages ====== (from eco925)')
  install.packages('dplyr')
}
if (!require('vroom')) {
  message('==== installing required packages ====== (from eco925)')
  install.packages('vroom')
}
if (!require('stringr')) {
  message('==== installing required packages ====== (from eco925)')
  install.packages('stringr')
}

library(dplyr)
library(vroom)

### cluster info
cluster_info = vroom(cluster)
colnames(cluster_info)=c('rep_seq','mem_seq')
# mem_seq = cluster_info %>% distinct(mem_seq) 
# rep_seq = cluster_info %>% distinct(rep_seq) 

sample_list = fs::dir_ls(input,type='file')
sample_name = strsplit(sample_list,'/') %>% rapply(., function(x) (tail(x,1))) # remove base dir
sample_name = strsplit(sample_name,'\\.') %>% rapply(., function(x) (head(x,1))) # remove tail

cluster_table = cluster_info
message('=========== start to merge table ==============')
for (i in 1:length(sample_list)) {
  dt = vroom(sample_list[i])
  dt = dt[,-2]
  colnames(dt)=c('mem_seq',sample_name[i])
  cluster_table = full_join(cluster_table,dt,by='mem_seq')
}

cluster_table = cluster_table[,-2]
cluster_table = cluster_table %>% group_by(rep_seq) %>% summarise(across(everything(), sum))
message('=========== start to write output ==============')
vroom_write(cluster_table, file='cluster_table.tsv')
