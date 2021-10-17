### we used bbmap tools to examine the output of different assembly strategies. Mainly, individual assembly or co-assembly
############################### co-assembly examination #####################################
##### bbmap tools #####
### take a try in files of L1
### estimate coverage after each steps
###
mkdir L1
ln -s ../L1-* .
ls -d *.R1.fq.gz | cut -d '.' -f1 > L1_list.txt
### bbmap
# Using BBMap, reads that were mapped to masked human, cat, dog, and mouse references at 95% identity and aligned to common microbial contaminants were separated.
# To index and map at the same time:
# To build an index in-memory without writing to disk:
# bbmap.sh in=${i}.qc.fq out=${i}.mapped.sam ref=ref.fa nodisk -t=40
# To split input into mapped and unmapped, in fastq format:
# bbmap.sh in=${i}.mapped.sam outu=u${i}.rmhos.fq -t=40

### bbduk
# used to remove contaminants, reads that contained adapter sequences, and right quality trim reads where quality dropped to 0.
# BBDuk was also applied to eliminate reads containing 1 or more “N” bases, having an average quality score across the read of less than 13 or containing a minimum length of ≤ 41 bp or 33% of the full read length.
# parameters used in reference: right end, Kmer length = 27, maximum substitutions = 1, minimum quality = 20, minimum overlap = 20, minimum length = 20
### bbduk
# adapter remove + qc control
parallel -j 10 --xapply 'bbduk.sh in1={1}.R1.fq.gz in2={1}.R2.fq.gz out={1}.qc.fq ref=~/bbmap_resources/adapters.fa ktrim=r k=23 mink=11 hdist=1  qtrim=r trimq=20 maq=20 -t=4 minlength=50 tpe tbo' :::: L1_list.txt
# remove phix protein (no-effect in metagenomics?)
parallel -j 11 --xapply 'bbduk.sh in={1}.qc.fq out={1}.clean.fq ref=~/bbmap_resources/phix.fa k=31 hdist=1 stats={1}.stats.txt -t=6' :::: L1_list.txt

### bbmerge
# It is designed for kmer-based operations using Tadpole, which include both merging overlapping and non-overlapping reads, kmer-based error-correction, and kmer-based filtering.
# Kmer-based operations should only be used with shotgun (randomly-fragmented) libraries, never with amplicon libraries (such as 16S).
#  If you run BBMerge, and under, say, 15% of the reads merge, even at very loose stringency, it’s probably a waste of time to merge – you’ll just make the workflow more complicated, and possibly get a lot of false-positives. Also, don’t try to merge single-ended libraries or long-mate-pair libraries that are not in an “innie” orientation.
# effect for assembly in spades by author:http://seqanswers.com/forums/showthread.php?t=43906&page=5
# Result: BBmerge looks unnecessary due to the only 20% merge rate and variated insert size of library (400-600bp in average)
# historical codes:
while read prefix; do
  bbmerge-auto.sh in=${prefix}.clean.fq out=${prefix}.bmer.fq outu=${prefix}.bumer.fq ecct extend2=20 iterations=5 -Xmx250g -t=40
done <L1_list.txt
# PS: see our binning log files.

### bbnorm
# bbnorm to adjust coverage will be good for spades:
# http://seqanswers.com/forums/showthread.php?t=49763
# loglog.sh in1=merged.R1.fq in2=merged.R2.fq # ~60G 5000s
bbnorm.sh in1=merged.R1.fq in2=merged.R2.fq out=highpass.fq outt=lowpass.fq passes=1 target=40 min=5 -Xmx250g threads=40 # consume 230G memory actually, ~3 hours

cd ..
mkdir 02_spades
cd 02_spades/L1
ln -s ../../01_bbmap/L1/*.clean.fq .
ln -s ../../01_bbmap/L1/L1_list.txt .
# parallel bbnorm on individual samples
parallel -j 2 --xapply 'bbnorm.sh in={}.clean.fq out={}.bnorm.fq passes=1 target=999999999 min=2 -Xmx120g threads=10' :::: L1_list.txt
# or
while read prefix; do
  bbnorm.sh in=${prefix}.clean.fq out=${prefix}.norm.fq passes=1 target=999999999 min=2 -Xmx120g -t=30
done <L1_list.txt
# for two-step which set to min =5
bbnorm.sh in=L1-all.tmp.fq out=test.L1-all.2step.fq passes=2 target=999999999 min=5 -Xmx200g -t=50
### we did not show detailed codes of normalizing co-assembly fastq, because it is similar with above

### Nonpareil coverage estimation
# nonpareil should start after qc trimming (> q20)
# https://nonpareil.readthedocs.io/en/latest/redundancy.html
# take a try
mkdir cov_est
# create qualified seq files for following test
seqkit seq -j 10 -m 50 L1-H1.clean.fq > test.L1-H1.clean.fq
seqkit seq -j 10 -m 50 L1-H1.bumer.fq > test.L1-H1.bumer.fq
cat L1*.bnorm.fq > L1-all.tmp.fq
bbnorm.sh in=L1-all.tmp.fq out=test.L1-all.2step.fq passes=2 target=999999999 min=5 -Xmx200g -t=50
### in following steps: co-assembly means we concatenated seq files together
# test individual fastq coverage (bbduk)
nonpareil -s test.L1-H1.clean.fq -T kmer -f fastq -b cov_est/clean -t 40 -R 250000
# test individual merge coverage (bbduk + bbmerge)
nonpareil -s L1-H1.bmer.fq -T kmer -f fastq -b cov_est/bmer -t 40 -R 250000
# test individual non-merge coverage (bbduk + bbmerge)
nonpareil -s test.L1-H1.bumer.fq -T kmer -f fastq -b cov_est/bumer -t 40 -R 100000
# test individual normalized coverage (bbduk + bbnorm)
seqkit seq -j 10 -m 50 L1-H1.bnorm.fq > test.L1-H1.bnorm.fq
nonpareil -s test.L1-H1.bnorm.fq -T kmer -f fastq -b cov_est/bnorm -t 40 -R 100000
# test co-assembly normalized coverage (bbduk + co-assembly + bbnorm)
nonpareil -s test.L1-all.bnorm.fq -T kmer -f fastq -b cov_est/co-bnorm -t 40 -R 300000
# test two-step normalized coverage (bbduk + bbnorm + co-assembly + bbnorm)


# PS: see our binning log files.
