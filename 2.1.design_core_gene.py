#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

from utils import (create_dir, find_closest_region, find_overlap, link_region,
                   read_core_gene_trans, read_noney)

# 对核心基因设计探针

# 核心基因转录本坐标文件
core_gene_region_file = 'core_gene_region.txt'
# 大区文件
noney_file = 'noneY.txt'
# 输出目录
noney_dir = "./noney"
core_gene_dir = os.path.join(noney_dir, 'core_gene')
output_file = os.path.join(core_gene_dir, "core_gene_design.txt")
create_dir(noney_dir)
create_dir(core_gene_dir)
# 起始、终止密码子上下扩展长度
expand = 1000

# 读入数据
noney = read_noney(noney_file)  # 数据存储格式 noney[bigzone]['row' + str(row)] = tmps
core_gene_region = read_core_gene_trans(core_gene_region_file)

def expand_n(name, chrom, center, expand):
    start = center - expand
    end = center + expand
    return name + "," + chrom + ":" + str(start) + "-" + str(end)


# 设计
designed_region = {}
for bigzone in noney.keys():
    # 核心基因
    genes = []
    # 核心区域
    critical_regions_raw = []
    critical_regions_raw_str = []
    for row in noney[bigzone].keys():
        genes.extend(noney[bigzone][row]['核心区关键基因'])
        if noney[bigzone][row]['critical_start'] != 0:
            chrom = noney[bigzone][row]['Chr']
            critical_start = noney[bigzone][row]['critical_start']
            critical_end = noney[bigzone][row]['critical_end']
            critical_regions_raw.append({"chrom": chrom, "start": critical_start, "end": critical_end })
            critical_regions_raw_str.append(chrom + ":" + str(critical_start) + "-" + str(critical_end))
    critical_regions = link_region(critical_regions_raw)
    # 去冗余
    genes = set(genes)
    # 设计
    for gene in genes:
        # 探针命名的前缀
        prefix = bigzone + "-" + gene
        if gene not in core_gene_region:
            print(f'[Error] {gene} 缺失转录本区域信息')
            continue
        chrom = core_gene_region[gene]['start_codon'][0]['chrom']
        strand = core_gene_region[gene]['start_codon'][0]['strand']
        # 1. 与核心区域取交
        # 起始密码子与核心区取交集 [{}, {}]
        start_codon_overlaps = find_overlap(core_gene_region[gene]['start_codon'][0], [critical_regions])
        # 终止密码子与核心区取交集 [{}, {}]
        stop_codon_overlaps = find_overlap(core_gene_region[gene]['stop_codon'][0], [critical_regions])
        # 编码区 覆盖区域 {}
        codon_region = link_region([core_gene_region[gene]['stop_codon'][0], core_gene_region[gene]['start_codon'][0]])
        # 编码区域与核心区取交集 [{}, {}]
        codon_overlap = find_overlap(codon_region, [critical_regions])
        # 编码区域与核心区的交集 连接成一个大区域
        codon_overlap_link = link_region(codon_overlap)
        # 找出交集区域与CDS每一个区域的重叠区域
        cds_overlap = find_overlap(codon_overlap_link, core_gene_region[gene]['CDS'])
        
        
        # 状态判断
        # 1. 没有核心区域 或者 与核心区没有任何交集
        status = '.'
        if critical_regions['chrom'] == "":
            # 没有核心区域：方案4设计
            status = 'no_critical'
        elif len(start_codon_overlaps) > 0  and len(stop_codon_overlaps) > 0:
            # 起始、终止都覆盖到了：方案1设计
            status = 'two_cover'
        elif len(start_codon_overlaps) == 0  and len(stop_codon_overlaps) > 0:
            # 只有终止密码子覆盖上了
            status = 'stop_cover'
        elif len(start_codon_overlaps) > 0  and len(stop_codon_overlaps) == 0:
            # 只有起始密码子覆盖上了
            status = 'start_cover'
        else:
            # 起始、终止密码子都没有覆盖到
            # 如果有交集: 方案3 设计
            if len(codon_overlap) > 0:
                status = 'inner_cover'
            else:
                # 如果没有交集： 方案4设计
                status = 'no_cover'

        
        # 根据状态，设计探针坐标区域
        # print(gene, status)
        if status == "no_critical" or status == 'two_cover' or status == 'no_cover':
            # 方式1：起始密码子、终止密码子、中心点exon各设计一个探针
            start = core_gene_region[gene]['start_codon'][0]['start']
            end = core_gene_region[gene]['stop_codon'][0]['end']
            if strand == '-':
                start = core_gene_region[gene]['start_codon'][0]['end']
                end = core_gene_region[gene]['stop_codon'][0]['start']
            point = int((start + end) / 2)
            closest_cds_point, closest_cds_point_candidate = find_closest_region({"chrom": chrom, "start":point, "end": point}, core_gene_region[gene]['CDS'], 2)
            
            name = prefix + "-1"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'start_codon',
                'middle_point': start,
                "chrom": chrom,
                "start": start  - expand,
                "end": start + expand,
                "candidate_exon": '.',
                "candidate_updown": expand_n(name + "-expand2K", chrom, start, 2000) + ";" + expand_n(name + "-expand5K", chrom, start, 5000)
            }
            name = prefix + "-2"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'closest_CDS',
                'middle_point': point,
                "chrom": chrom,
                "start": closest_cds_point['start'],
                "end": closest_cds_point['end'],
                "candidate_exon": closest_cds_point_candidate,
                "candidate_updown": "."
            }
            name = prefix + "-3"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'stop_codon',
                'middle_point': end,
                "chrom": chrom,
                "start": end - expand,
                "end": end + expand,
                "candidate_exon": ".",
                "candidate_updown": expand_n(name + "-expand2K", chrom, end, 2000) + ";" + expand_n(name + "-expand5K", chrom, end, 5000)
            }
        elif status == "start_cover":
            # 方式2：起始密码子覆盖
            start = core_gene_region[gene]['start_codon'][0]['start']
            demarcation_point = codon_overlap_link['end'] # 分界点
            if strand == '-':
                start = core_gene_region[gene]['start_codon'][0]['end']
                demarcation_point = codon_overlap_link['start']
            # 中心点
            point = int((start + demarcation_point) / 2)
            # 1. 分界点最近的cds
            closest_cds_demarcation_point, closest_cds_demarcation_point_candidate = find_closest_region({"chrom": chrom, "start":demarcation_point, "end": demarcation_point}, cds_overlap, 2)
            # 2. 与中心点最近的cds
            closest_cds_point, closest_cds_point_candidate = find_closest_region({"chrom": chrom, "start":point, "end": point}, cds_overlap, 2)

            name = prefix + "-1"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'start_codon',
                'middle_point': start,
                "chrom": chrom,
                "start": start  - expand,
                "end": start + expand,
                "candidate_exon": ".",
                "candidate_updown": expand_n(name + "-expand2K", chrom, start, 2000) + ";" + expand_n(name + "-expand5K", chrom, start, 5000)
            }
            name = prefix + "-2"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'closest_CDS',
                'middle_point': point,
                "chrom": chrom,
                "start": closest_cds_point['start'],
                "end": closest_cds_point['end'],
                "candidate_exon": closest_cds_point_candidate,
                "candidate_updown": "."
            }
            name = prefix + "-3"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'demarcation_point',
                'middle_point': demarcation_point,
                "chrom": chrom,
                "start": closest_cds_demarcation_point['start'],
                "end": closest_cds_demarcation_point['end'],
                "candidate_exon": closest_cds_demarcation_point_candidate,
                "candidate_updown": "."
            }
        elif status == "stop_cover":
            # 方式3：终止密码子覆盖
            end = core_gene_region[gene]['stop_codon'][0]['end']
            demarcation_point = codon_overlap_link['start'] # 分界点
            if strand == '-':
                end = core_gene_region[gene]['stop_codon'][0]['start']
                demarcation_point = codon_overlap_link['end']
            # 中心点
            point = int((end + demarcation_point) / 2)
            # 1. 分界点最近的cds
            closest_cds_demarcation_point, closest_cds_demarcation_point_candidate = find_closest_region({"chrom": chrom, "start":demarcation_point, "end": demarcation_point}, cds_overlap, 2)
            # 2. 与中心点最近的cds
            closest_cds_point, closest_cds_point_candidate = find_closest_region({"chrom": chrom, "start":point, "end": point}, cds_overlap, 2)

            name = prefix + "-1"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'demarcation_point',
                'middle_point': demarcation_point,
                "chrom": chrom,
                "start": closest_cds_demarcation_point['start'],
                "end": closest_cds_demarcation_point['end'],
                "candidate_exon": closest_cds_demarcation_point_candidate,
                "candidate_updown": "."
            }
            name = prefix + "-2"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'closest_CDS',
                'middle_point': point,
                "chrom": chrom,
                "start": closest_cds_point['start'],
                "end": closest_cds_point['end'],
                "candidate_exon": closest_cds_point_candidate,
                "candidate_updown": "."
            }
            name = prefix + "-3"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'stop_codon',
                'middle_point': end,
                "chrom": chrom,
                "start": end - expand,
                "end": end + expand,
                "candidate_exon": ".",
                "candidate_updown": expand_n(name + "-expand2K", chrom, end, 2000) + ";" + expand_n(name + "-expand5K", chrom, end, 5000)
            }
        elif status == "inner_cover":
            # 方式3： 内部被包含
            demarcation_point_start = codon_overlap_link['start']
            demarcation_point_end = codon_overlap_link['end']
            point = int((demarcation_point_start + demarcation_point_end) / 2)
            if strand == '-':
                demarcation_point_start, demarcation_point_end = demarcation_point_start, demarcation_point_end
            # 1. start cds
            closest_cds_start, closest_cds_start_candidate = find_closest_region({"chrom": chrom, "start":demarcation_point_start, "end": demarcation_point_start}, cds_overlap, 2)
            # 2. center cds
            closest_cds_point, closest_cds_point_candidate = find_closest_region({"chrom": chrom, "start":point, "end": point}, cds_overlap, 2)
            # 3. end cds
            closest_cds_end, closest_cds_end_candidate = find_closest_region({"chrom": chrom, "start":demarcation_point_end, "end": demarcation_point_end}, cds_overlap, 2)
            
            name = prefix + "-1"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'demarcation_point',
                'middle_point': demarcation_point_start,
                "chrom": chrom,
                "start": closest_cds_start['start'],
                "end": closest_cds_start['end'],
                "candidate_exon": closest_cds_start_candidate,
                "candidate_updown": "."
            }
            name = prefix + "-2"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'closest_CDS',
                'middle_point': point,
                "chrom": chrom,
                "start": closest_cds_point['start'],
                "end": closest_cds_point['end'],
                "candidate_exon": closest_cds_point_candidate,
                "candidate_updown": "."
            }
            name = prefix + "-3"
            designed_region[name] = {
                "name": name,
                'status': status,
                "type": 'demarcation_point',
                'middle_point': demarcation_point_end,
                "chrom": chrom,
                "start": closest_cds_end['start'],
                "end": closest_cds_end['end'],
                "candidate_exon": closest_cds_end_candidate,
                "candidate_updown": "."
            }
        
        # 追加额外注释
        for name in [prefix + "-1", prefix + "-2", prefix + "-3"]:
            designed_region[name]['gene'] = gene
            designed_region[name]['mrna'] = core_gene_region[gene]['mrna']
            designed_region[name]['cds_count'] = len(core_gene_region[gene]['CDS'])
            designed_region[name]['strand'] = strand
            designed_region[name]['width'] = designed_region[name]['end'] - designed_region[name]['start'] + 1
            designed_region[name]['center'] = int((designed_region[name]['end'] + designed_region[name]['start']) / 2)
            designed_region[name]['cnvplex_format'] = name + "," + designed_region[name]['chrom'] + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
            if designed_region[name]['candidate_exon'] != '.':
                designed_region[name]['candidate_exon'] =  ";".join([ name + "-candidate_" + exon for exon in designed_region[name]['candidate_exon'].split(";") ])
            designed_region[name]['critical_region_raw'] = ",".join(critical_regions_raw_str) if critical_regions['chrom'] != "" else "."
            designed_region[name]['critical_region_merge'] = critical_regions['chrom'] + ":" + str(critical_regions['start']) + "-" + str(critical_regions['end']) if critical_regions['chrom'] != "" else "."
            designed_region[name]['cds'] = ",".join([ region['chrom'] + ":" + str(region['start']) + "-" + str(region['end']) for region in core_gene_region[gene]['CDS']])
            designed_region[name]['start_codon'] =  core_gene_region[gene]['start_codon'][0]['chrom'] + ":" + str(core_gene_region[gene]['start_codon'][0]['start']) + "-" + str(core_gene_region[gene]['start_codon'][0]['end'])
            designed_region[name]['stop_codon'] = core_gene_region[gene]['stop_codon'][0]['chrom'] + ":" + str(core_gene_region[gene]['stop_codon'][0]['start']) + "-" + str(core_gene_region[gene]['stop_codon'][0]['end'])
        
# 输出
titles = ['name', 'gene', 'mrna', 'cds_count', 'strand', 'status', 'type', 'middle_point', 'center', 'chrom', 'start', 'end', 'width', 'candidate_exon', 'candidate_updown', 'cnvplex_format', 'critical_region_raw', 'critical_region_merge', 'cds', 'start_codon', 'stop_codon']
with open(output_file, 'w') as fh:
    fh.write("\t".join(titles) + "\n")
    for name in designed_region.keys():
        values = [str(designed_region[name][title]) for title in titles]
        fh.write("\t".join(values) + "\n")
print("输出文件：" + output_file)