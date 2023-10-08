#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

import pybedtools

# 找到探针覆盖范围内的unsoftmask区域，目的：降低宽度。减少cnvplex运行时间消耗
input_dir = "./noney/sub_region/"
input_file = os.path.join(input_dir, "sub_region_design.txt") # 探针设计文件
unsoftmask_file = './target.txt.nocnv'  # 非小写字母覆盖区域
output_file = os.path.join(input_dir, "sub_region_design.unsoftmask.txt")


def read_input(file):
    data = {}
    with open(file, 'r', encoding='unicode_escape') as fh:
        heads = fh.readline().rstrip().split('\t')
        for line in fh:
            values = line.strip().split('\t')
            tmps = { heads[col]: values[col]  for col in range(0,len(values))}
            data[tmps['name']] = tmps
    return data


# 读入数据
input_data = read_input(input_file)
unsoftmask = pybedtools.BedTool(unsoftmask_file)

# 每个步长的数据计算overlap
for step in ['cnvplex_1K', 'cnvplex_10K', 'cnvplex_50K', 'cnvplex_100K', 'cnvplex_200K']:
    print("处理步长：" + str(step))
    
    # 生成步长bed文件
    step_file = os.path.join(input_dir, step + ".bed")
    with open(step_file, 'w') as fh:
        for name in input_data.keys():
            if input_data[name][step] == '.':
                continue
            chrom, start, end = re.split('[:-]', input_data[name][step].split(',')[1]) 
            fh.write("\t".join([chrom, start, end, name]) + "\n")

    # 步长区域内，找 unsoftmask region
    step_obj = pybedtools.BedTool(step_file)
    finds = {}
    for field in step_obj.intersect(unsoftmask, wao=True):
        name = field.fields[3]
        p_chrom, p_start, p_end, chrom, start, end, tlen = field.fields[0], field.fields[1], field.fields[2], field.fields[4], field.fields[5], field.fields[6], field.fields[8]
        # 重叠区域至少10bp
        if tlen == "." or tlen == '0':
            continue
        else:
            tlen = int(tlen)
        if tlen < 10:
            continue
        # 做记录
        if name not in finds:
            finds[name] = {'target': [], 'length': 0, 'count': 0}
        # 找交集区域
        start_overlap = max(int(p_start), int(start))
        end_overlap = min(int(p_end), int(end))
        
        finds[name]['target'].append(name + "," + chrom + ":" + str(start_overlap) + "-" + str(end_overlap))
        finds[name]['length'] += tlen
        finds[name]['count'] += 1
    
    # 放入数据库
    for name in input_data.keys():
        if name not in finds:
            input_data[name][step + "_unsoftmask"] = '.'
            input_data[name][step + "_unsoftmask_count"] = '.'
            input_data[name][step + "_unsoftmask_length"] = '.'
        else:
            input_data[name][step + "_unsoftmask"] = ";".join(finds[name]['target'])
            input_data[name][step + "_unsoftmask_count"] = finds[name]['count']
            input_data[name][step + "_unsoftmask_length"] = finds[name]['length']


# 输出
titles = []
with open(input_file, 'r') as fh:
    titles = fh.readline().strip().split('\t')

for step in ['cnvplex_1K', 'cnvplex_10K', 'cnvplex_50K', 'cnvplex_100K', 'cnvplex_200K']:
    titles.append(step + "_unsoftmask")
    titles.append(step + "_unsoftmask_count")
    titles.append(step + "_unsoftmask_length")

with open(output_file, 'w') as fh:
    fh.write("\t".join(titles) + "\n")
    for name in input_data.keys():
        values = [str(input_data[name][title]) for title in titles]
        fh.write("\t".join(values) + "\n")
print("输出文件：" + output_file)