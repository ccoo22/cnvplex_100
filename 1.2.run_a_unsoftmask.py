#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

import pybedtools

# 找到a类探针覆盖范围内的unsoftmask区域，目的：降低宽度。减少cnvplex运行时间消耗
input_file = 'chromosome/final.txt'  # a类探针设计文件
unsoftmask_file = './target.txt'  # 非小写字母覆盖区域
input_dir = "./chromosome"
output_file = "chromosome/final.txt.unsoftmask.txt"
steps = [10, 50, 100, 200] # 每次跨越步长


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
for step in steps:
    need_title = 'expand_' + str(step) + "K"
    print("处理步长：" + str(step))
    
    # 生成步长bed文件
    step_file = os.path.join(input_dir, "step"+str(step) + ".bed")
    with open(step_file, 'w') as fh:
        for name in input_data.keys():
            chrom, start, end = re.split('[:-]', input_data[name][need_title].split(',')[1]) 
            fh.write("\t".join([chrom, start, end, name]) + "\n")

    # 步长区域内，找 unsoftmask region
    step_obj = pybedtools.BedTool(step_file)
    finds = {}
    for field in step_obj.intersect(unsoftmask, wao=True):
        name = field.fields[3]
        chrom, start, end, tlen = field.fields[4], field.fields[5], field.fields[6], field.fields[7]
        # 重叠区域至少100bp
        if tlen == "." or tlen == '0':
            continue
        else:
            tlen = int(tlen)
        if tlen < 100:
            continue
        # 做记录
        if name not in finds:
            finds[name] = {'target': [], 'length': 0, 'count': 0}
        finds[name]['target'].append(name + "," + chrom + ":" + start + "-" + end)
        finds[name]['length'] += tlen
        finds[name]['count'] += 1
    
    # 放入数据库
    for name in input_data.keys():
        if name not in finds:
            input_data[name][need_title + "_unsoftmask"] = '.'
            input_data[name][need_title + "_unsoftmask_count"] = '.'
            input_data[name][need_title + "_unsoftmask_length"] = '.'
        else:
            input_data[name][need_title + "_unsoftmask"] = ";".join(finds[name]['target'])
            input_data[name][need_title + "_unsoftmask_count"] = finds[name]['count']
            input_data[name][need_title + "_unsoftmask_length"] = finds[name]['length']


# 输出
titles = ['name', 'chrom', 'center', 'type']
for step in steps:
    titles.append('expand_' + str(step) + "K")
    titles.append('expand_' + str(step) + "K" + "_unsoftmask")
    titles.append('expand_' + str(step) + "K" + "_unsoftmask_count")
    titles.append('expand_' + str(step) + "K" + "_unsoftmask_length")
with open(output_file, 'w') as fh:
    fh.write("\t".join(titles) + "\n")
    for name in input_data.keys():
        values = [str(input_data[name][title]) for title in titles]
        fh.write("\t".join(values) + "\n")
print("输出文件：" + output_file)