#！/usr/bin/python3.6
# -*- coding:utf-8 -*-
# use like: tag_split merge.fastq split_dir

import re
import os
import sys

def mkdir(path): # waiting for be moved into base tools
	folder = os.path.exists(path)
	if not folder:
		os.makedirs(path)
		print ("--- Create a new directory ---")
	else:
		print ("--- The directory has already existed! ---")

# saved path (dir)
fdir = sys.argv[2]
mkdir(fdir)

# re match tags
tag = input('Please tell me the tag of sequences: ')

# I/O in low-memory consuming
with open('{}'.format(sys.argv[1]),'r') as file:
    for line in file:
        sample = line[line.index(tag):].strip(tag).strip()
        fastq = open('{}/{}.fastq'.format(fdir,sample),'a')
        fastq.write(line[:line.index(tag)]+'\n')
        for i in range(3): fastq.write(next(file))
file.close()
