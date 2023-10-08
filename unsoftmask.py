#!/home/genesky/software/python/3.9.4/bin/python3
import argparse
import logging
import os
import re
import sys

import pybedtools


def set_and_parse_args():
    """参数解析"""
    parser = argparse.ArgumentParser(
        description="获取cnvplex数据的unsoftmask区域", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--input', '-i', type=str, required=True,
                        help="输入文件，一列，每一行的数据格式为 name,chrom:start-end")
    args = parser.parse_args()

    return args


args = set_and_parse_args()


# 找到探针覆盖范围内的unsoftmask区域，目的：降低宽度。减少cnvplex运行时间消耗
input_file = args.input
unsoftmask_file = './target.txt.nocnv'
step_file = input_file + ".bed"
output_file = input_file + ".unsoftmask.txt"

def read_input(file):
    data = []
    with open(file, 'r', encoding='unicode_escape') as fh:
        for line in fh:
            if re.search('\w', line):
                data.append(line.strip())
    return data


# 读入数据
input_data = read_input(input_file)
unsoftmask = pybedtools.BedTool(unsoftmask_file)


# 生成步长bed文件

with open(step_file, 'w') as fh:
    for data in input_data:
        name, region = data.split(',')
        chrom, start, end = re.split('[:-]', region) 
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
    
    finds[name]['target'].append(name + "," + chrom + ":" + str(start_overlap) + "-" + str(end_overlap) + '\t' + str(end_overlap - start_overlap))
    finds[name]['length'] += tlen
    finds[name]['count'] += 1

# 输出

with open(output_file, 'w') as fh:
    for name in finds.keys():
        fh.write("\n".join(finds[name]['target']) + "\n")

print("结果：" + output_file)
