import numpy as np
import pandas as pd
import os,time,argparse,sys
#########################
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
This scripts is to merge coverage table of bamm by cluster info \n
please placed covergae file in their respective folder. \n
2021-09\n
""",formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i','--input', required=True,help='input folder')
    parser.add_argument('-c','--cluster', required=True,help='tsv folder')
    parser.add_argument('-o','--output',default='output.tsv',help='output merge table')
    return parser

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    input = args.input
    cluster = args.cluster
    output = args.output

sample_list = os.listdir(input)
dt=pd.read_csv(cluster,sep='\t',header=None,names=['Rep_seq','Seq_ID'])

print('------{}: start to process coverage tables'.format(time.strftime("%Y-%m-%d %X")))
for cov_path in sample_list:
    sample_name  = cov_path[:cov_path.index('.')]
    cov_path = os.path.join(input,cov_path)
    cov_dt = pd.read_csv(cov_path,sep='\t',header=0, usecols=[0,2],names=["Seq_ID", sample_name])
    dt = pd.merge(dt,cov_dt,on='Seq_ID',how='left')

print('------{}: data is merging'.format(time.strftime("%Y-%m-%d %X")))
index_dt = dt.set_index(["Seq_ID"])

merge_dt = grouped = index_dt.groupby('Rep_seq').agg('sum')

merge_dt.to_csv(output, sep='\t')
print('------{}: done'.format(time.strftime("%Y-%m-%d %X")))