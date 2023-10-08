#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

from utils import (create_dir, design_75_1, design_75_2, design_cnvplex,
                   read_noney, region_probe_design)

# 对核心区域设计探针

# 大区文件
noney_file = 'noneY.txt'
# 输出目录
noney_dir = "./noney"
core_gene_dir = os.path.join(noney_dir, 'sub_region')
output_file = os.path.join(core_gene_dir, "sub_region_design.txt")
create_dir(noney_dir)
create_dir(core_gene_dir)


# 读入数据
noney = read_noney(noney_file)  # 数据存储格式 noney[bigzone]['row' + str(row)] = tmps

extra_titles = ["L-REGION50", "L-REGION75", "R-REGION75", "R-REGION50", "REGION90 Sample Number"]

# 设计
designed_region = {}
for bigzone in noney.keys():
    for row in noney[bigzone].keys():
        # 对每一个亚区进行设计
        # 没有亚区的，不用设计
        if noney[bigzone][row]['REGION90-P1'] == 0:
            continue
        # 1. 亚区
        chrom = noney[bigzone][row]['Chr']
        start = noney[bigzone][row]['REGION90-P1']
        end = noney[bigzone][row]['REGION90-P2']
        width = end - start + 1
        target_regions = region_probe_design(chrom, start, end, 0.3)  # surfix  interval chrom center start end width probe_number
        for target_region in target_regions:
            name = bigzone + "-" + chrom + "_" + str(start) + "_" + str(end) + "-" + target_region['surfix']
            designed_region[name] = {}
            designed_region[name]['name'] = name
            designed_region[name]['sub_chrom'] = chrom 
            designed_region[name]['sub_start'] = start
            designed_region[name]['sub_end'] = end
            designed_region[name]['sub_width'] = end - start + 1
            designed_region[name]['probe_number'] = target_region['probe_number']
            designed_region[name]['interval'] = target_region['interval']
            designed_region[name]['center'] = target_region['center']
            designed_region[name]['chrom'] = target_region['chrom']
            designed_region[name]['start'] = target_region['start']
            designed_region[name]['end'] = target_region['end']
            designed_region[name]['width'] = target_region['width']
            designed_region[name]['exclude_if'] = target_region['exclude_if']
            # 补充额外title
            for extra_title in extra_titles:
                designed_region[name][extra_title] = noney[bigzone][row][extra_title]

        # 2. 左右两侧 50%、75% 区域
        side_names = []
        for side in ['L', 'R']:  # 分左右
            start_50 = noney[bigzone][row][f'{side}-REGION50']
            start_75 = noney[bigzone][row][f'{side}-REGION75']
            if start_50 != 0:
                info = {"75": {}, "50": {}}
                width_75 = start - start_75 if side == "L" else start_75 - end
                width_50 = start_75 - start_50 if side == "L" else start_50 - start_75
                # 是否设计探针？
                if width_75 > width * 0.1:
                    info['75'] = {"chrom": chrom, "start": start_75, 'width': width_75, "status": "75"}
                    if width_50 > (width + width_75) * 0.1:
                        info['50'] = {"chrom": chrom, "start": start_50, 'width': width_50, "status": '50'}
                else:
                    if (width_75 + width_50) > width * 0.1:
                        info['50'] = {"chrom": chrom, "start": start_50, 'width': width_75 + width_50, "status": "50_75"}
                # 根据 75% 50% 探针的设计原则，设计探针位置
                if info['75']:
                    if info['75']['width'] > width * 2:
                        # 大于核心区域2倍，则设计2个探针
                        prefix = bigzone + "-" + chrom + "_" + str(start) + "_" + str(end) + "-" + side + info['75']['status']
                        side_names.append(prefix + "-1")
                        side_names.append(prefix + "-2")
                        design_75_2(prefix, designed_region, side, info['75'])
                    else:
                        # 只需1个探针
                        name = bigzone + "-" + chrom + "_" + str(start) + "_" + str(end) + "-" + side + info['75']['status']
                        side_names.append(name)
                        design_75_1(name, designed_region, side, info['75'])
                if info['50']:
                    name = bigzone + "-" + chrom + "_" + str(start) + "_" + str(end) + "-" + side + info['50']['status']
                    side_names.append(name)
                    design_75_1(name, designed_region, side, info['50'])  # 它的设计规则和 75_1 完全一样
        # 补充其他信息
        for side_name in side_names:
            designed_region[side_name]['name'] = side_name
            designed_region[side_name]['sub_chrom'] = chrom 
            designed_region[side_name]['sub_start'] = start
            designed_region[side_name]['sub_end'] = end
            designed_region[side_name]['sub_width'] = end - start + 1
            designed_region[side_name]['interval'] = '.'
            designed_region[side_name]['chrom'] = chrom
            # 补充额外title
            for extra_title in extra_titles:
                designed_region[side_name][extra_title] = noney[bigzone][row][extra_title]

# 设计cnvplex格式区域，方便使用
design_cnvplex(designed_region)

# 输出
titles = ['name', 'sub_chrom', 'sub_start', 'sub_end', 'sub_width']
titles.extend(extra_titles)
titles.extend(['probe_number', 'interval', 'center', "chrom", 'start', 'end', 'width', 'exclude_if', 'cnvplex_1K', 'cnvplex_10K', 'cnvplex_50K', 'cnvplex_100K', 'cnvplex_200K'])
with open(output_file, 'w') as fh:
    fh.write("\t".join(titles) + "\n")
    for name in designed_region.keys():
        values = [str(designed_region[name][title]) for title in titles]
        fh.write("\t".join(values) + "\n")
print("输出文件：" + output_file)