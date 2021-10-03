#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Yufei,2021
# github.com/Yflyer

import time
import os,re,argparse,sys
import pandas as pd
table_prokka=sys.argv[1]
table_KEGG=sys.argv[2]
output=sys.argv[3]

print('------{}: data is loading'.format(time.strftime("%Y-%m-%d %X")))
dt1 = pd.read_csv(table_prokka, sep='\t', header=0)
dt2 = pd.read_csv(table_KEGG, sep='\t', header=0,low_memory=False)
dt2 = dt2.drop(index=0,columns='#')
print('------{}: data is merging'.format(time.strftime("%Y-%m-%d %X")))
df = pd.merge(dt1,dt2,left_on='Locus_ID', right_on='gene name',how='outer')
df
df.to_csv(output, sep='\t')
print('------{}: done'.format(time.strftime("%Y-%m-%d %X")))