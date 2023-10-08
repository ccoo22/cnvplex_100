#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

from utils import create_dir, design_cnvplex, read_noney, region_probe_design

# 对核心区域设计探针

# 大区文件
noney_file = 'noneY.txt'
# 输出目录
noney_dir = "./noney"
critical_region_dir = os.path.join(noney_dir, 'critical_region')
output_file = os.path.join(critical_region_dir, "critical_region_design.txt")
create_dir(noney_dir)
create_dir(critical_region_dir)


# 读入数据
noney = read_noney(noney_file)  # 数据存储格式 noney[bigzone]['row' + str(row)] = tmps

# 设计
designed_region = {}
for bigzone in noney.keys():
    for row in noney[bigzone].keys():
        # 对大区的每一行的核心区域单独设计
        # 没有核心区域的，不用设计
        if noney[bigzone][row]['critical_start'] == 0:
            continue
        chrom = noney[bigzone][row]['Chr']
        start = noney[bigzone][row]['critical_start']
        end = noney[bigzone][row]['critical_end']
        target_regions = region_probe_design(chrom, start, end, 0.2)  # surfix  interval chrom center start end width probe_number, 设计范围：center上下20%
        for target_region in target_regions:
            name = bigzone + "-" + chrom + "_" + str(start) + "_" + str(end) + "-" + target_region['surfix']
            designed_region[name] = {}
            designed_region[name]['name'] = name
            designed_region[name]['critical_chrom'] = chrom 
            designed_region[name]['critical_start'] = start
            designed_region[name]['critical_end'] = end
            designed_region[name]['critical_width'] = end - start + 1
            designed_region[name]['probe_number'] = target_region['probe_number']
            designed_region[name]['interval'] = target_region['interval']
            designed_region[name]['center'] = target_region['center']  # 用于计算探针与中心点的距离
            designed_region[name]['chrom'] = target_region['chrom']
            designed_region[name]['start'] = target_region['start']
            designed_region[name]['end'] = target_region['end']
            designed_region[name]['width'] = target_region['width']
            designed_region[name]['exclude_if'] = target_region['exclude_if']

# 制作cnvplex 格式数据
design_cnvplex(designed_region)


# 输出
titles = ['name', 'critical_chrom', 'critical_start', 'critical_end', 'critical_width', 'probe_number', 'interval', 'center', "chrom", 'start', 'end', 'width', 'exclude_if', 'cnvplex_1K', 'cnvplex_10K', 'cnvplex_100K', 'cnvplex_200K']
with open(output_file, 'w') as fh:
    fh.write("\t".join(titles) + "\n")
    for name in designed_region.keys():
        values = [str(designed_region[name][title]) for title in titles]
        fh.write("\t".join(values) + "\n")
print("输出文件：" + output_file)