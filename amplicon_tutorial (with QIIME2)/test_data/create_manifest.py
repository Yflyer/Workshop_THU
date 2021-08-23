import re
import os
import sys

# use listdir to get files names
# e.g.,in 'X2016.1N',X is the tag. the length of sample names are not equal, so use a tag rather a fixed index are better.
file1 = sys.argv[1]
file2 = sys.argv[2]
tag =  sys.argv[3]
f1paths = os.listdir(file1)
f2paths = os.listdir(file2)

# use index get the sample names
samples = [elem[elem.index(tag):] for elem in f1paths]

manifest = open ('manifest.tsv','w')
manifest.write('sample-id'+'\t'+'forward-absolute-filepath'+'\t'+'reverse-absolute-filepath'+'\n')
for (sample,f1path,f2path) in zip(samples,f1paths,f2paths):
    manifest.write(sample+'\t'+'$PWD/'+file1+'/'+f1path+'\t'+'$PWD/'+file2+'/'+f2path+'\n')
