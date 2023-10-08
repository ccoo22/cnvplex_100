#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

from utils import parse_gtf, read_core_gene

# 获取核心基因的转录本区域结构

# 染色体名称转换
chrom_ncbi_to_ucsc = {}
with open('GRCh37_latest_assembly_report.txt', 'r') as fh:
    for line in fh:
        if re.match('#', line):
            continue
        values = line.strip().split('\t')
        chrom_ncbi_to_ucsc[values[6]] = values[9]

gtf_file = '/home/pub/output2/research_and_customized_project/cnvplex_113/GCF_000001405.25_GRCh37.p13_genomic.gtf'
# 核心基因文件
core_gene_file = './core_gene.txt'
core_gene_region = "./core_gene_region.txt"

core_gene = read_core_gene(core_gene_file)
gtf_info = parse_gtf(gtf_file)
features = ['exon', 'CDS', 'start_codon', 'stop_codon']

with open(core_gene_region, 'w') as fh:
    fh.write("\t".join(['gene', 'trans', 'chrom', 'start', 'end', 'strand', 'feature', 'feature_number']) + "\n")
    for gene in core_gene.keys():
        trans = core_gene[gene]
        if trans not in gtf_info:
            fh.write("\t".join([gene, trans, '.', '.', '.', '.', '.', '.']) + "\n")
            continue
        else:
            for feature in features:
                # gtf_info[mrna][feature].append({'chrom': chrom, 'start': start, 'end': end, 'feature': feature, 'feature_number': feature_number})
                if feature in gtf_info[trans]:
                    for region in gtf_info[trans][feature]:
                        fh.write("\t".join([gene, trans, chrom_ncbi_to_ucsc[region['chrom']], str(region['start']), str(region['end']), region['strand'] ,region['feature'], region['feature_number']]) + "\n")