#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2021
# github.com/Yflyer
import time
import os,re,argparse,sys
from io import StringIO
import pandas as pd

"""Astror's toolkit for Eco925
github.com/Yflyer/
----------
file_name : reformat_prokka.py
This scripts is to reformat the tbl and tsv files of prokka
please placed these file in their respective folder.
2021-09
"""

def get_parser():
    parser = argparse.ArgumentParser(description="""Astror's toolkit for Eco925\n
github.com/Yflyer/ \n
----------
file_name : reformat_prokka.py \n
This scripts is to reformat the tbl and tsv files of prokka \n
please placed these file in their respective folder. \n
2021-09\n
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-l','--tbl_dir', required=True,help='tbl folder')
    parser.add_argument('-t','--tsv_dir', required=True,help='tsv folder')
    parser.add_argument('-o','--output',default='output.tsv',help='reformat result')
    return parser

def get_tbl_dict(input_f):
    dict_tbl={}
    for line in input_f:
        element = line.strip()
        if element[0] == '>':
            contig_ID=element[element.index(' '):]
            dict_tbl[contig_ID]=''
        else:
            dict_tbl[contig_ID]= dict_tbl[contig_ID]+'\t'+element+'\t'

    return dict_tbl

def get_tbl_info(sample_name,dict_tbl):
    tbl_info=''
    try:
        for contig_ID in dict_tbl.keys():
            len_info = re.findall(r"\s(\d+)\s(\d+)\s(\w+)",dict_tbl[contig_ID]) # finall start and end info in length
            locus_list = re.findall(r"locus_tag\t(\w+)",dict_tbl[contig_ID]) # findall all locus_tag
            if len(locus_list)>0:
                for i in range(len(locus_list)):
                    info = sample_name+'\t'+contig_ID+'\t'+locus_list[i]+'\t'+len_info[i][0]+'\t'+len_info[i][1]+'\t'+len_info[i][2]+'\n'
                    #output_f.write(info)
                    tbl_info=tbl_info+info
    except Exception as e:
        print('---------',e,':sample {} tbl can not processed'.format(sample_name),'---------')
    return tbl_dt

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    tbl_dir = args.tbl_dir
    tsv_dir = args.tsv_dir
    output = args.output


tbl_path=os.path.join(os. getcwd(),tbl_dir)
tbl_list = os.listdir(tbl_path)

tsv_path=os.path.join(os. getcwd(),tsv_dir)
tsv_list = os.listdir(tsv_path)

print('------{}: start to process tbl and tsv files from prokka'.format(time.strftime("%Y-%m-%d %X")))
tbl_dt='Sample_name\tContig_ID\tLocus_ID\tStart\tEnd\tType\n'
### extract info from tbl by dict and them storage by stringIO
for tbl_file in tbl_list:
    print('------{}: now start tbl file {}'.format(time.strftime("%Y-%m-%d %X"),tbl_file))
    sample_name  = tbl_file[:tbl_file.index('.')]
    with open (os.path.join(tbl_path,tbl_file),'r') as input_f:
        dict_tbl = get_tbl_dict(input_f)
        tbl_dt = tbl_dt+get_tbl_info(sample_name,dict_tbl)

index_dt = pd.DataFrame(pd.read_csv(StringIO(tbl_dt), sep='\t', header=0))

anno_dt=pd.DataFrame()
for tsv_file in tsv_list:
    sample_name  = tsv_file[:tsv_file.index('.')]
    print('------{}: now start tsv file {}'.format(time.strftime("%Y-%m-%d %X"),tsv_file))
    sample_dt = pd.DataFrame(pd.read_csv (os.path.join(tsv_path,tsv_file), sep='\t', header=0))
    sample_dt['Sample_name']=sample_name
    anno_dt = pd.concat([anno_dt,sample_dt])

print('------{}: data is merging'.format(time.strftime("%Y-%m-%d %X")))
anno_dt['index']=anno_dt["Sample_name"] + anno_dt["locus_tag"] + anno_dt["ftype"]
index_dt['index']=index_dt["Sample_name"] + index_dt["Locus_ID"] + index_dt["Type"]
df = pd.merge(index_dt,anno_dt[['index','gene','length_bp','product','EC_number','COG']],on='index', how='left')
df.to_csv(output, sep='\t')
print('------{}: done'.format(time.strftime("%Y-%m-%d %X")))
