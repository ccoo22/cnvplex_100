#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

import pysam

# 获取基因组非softmask区域的bed
genome_fa = '/home/genesky/database_new/ucsc/fasta/hg19/samtools_index/hg19.fa'

fh = open("target.txt", 'w')
ob = pysam.FastxFile(genome_fa)
for r in ob:
    chrom = r.name
    for iter in re.finditer('[ATCG]+', r.sequence, flags=0):
        start, end = iter.span()
        len = end - start 
        # # 区间太小的，不要，很难设计引物
        # if len < 100:
        #     continue
        fh.write(chrom + "\t" + str(start + 1) + "\t" + str(end) + "\t" + str(len) + "\n")
ob.close()
fh.close()