#!/home/genesky/software/python/3.9.4/bin/python3
import glob
import logging
import math
import multiprocessing
import os
import re
import sys
import warnings

import xlsxwriter
from utils import fix_by_manual, skip_C_5M

warnings.filterwarnings('ignore')
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout,
    level=logging.INFO
    )
log = logging.getLogger(__name__)

use_all = False
if len(sys.argv) > 1:
    use_all = True

# 把目标信息 与 探针信息组合起来
target_file ='./chromosome/final.txt.unsoftmask.txt'
probe_file = './probe/chromosome/final.txt' if not use_all else './probe/chromosome/final.txt.all.txt'
output_file = "./probe/chromosome/A_probe.xlsx"
output_file_raw = "./probe/chromosome/A_probe_raw.xlsx"
min_dist = 200 * 1000 # 目标位点与探针之间的最远距离
bestTm = 65 # 最佳Tm 值
probe_titles = ['chrom', 'start', 'end', 'strand', '最高同源性', 'CNV%', 'SNP（MAF）', '5_seq', '5_length', '5_tm', '3_seq', '3_length', '3_tm']

# 后期人工修正
fix_file = "./buchong/chromsome.txt"

print("探针设计目标区域文件：" + target_file)
print("探针文件：" + probe_file)


def read_input(file, key_title):
    data = {}
    with open(file, 'r') as fh:
        heads = fh.readline().rstrip().split('\t')

        for line in fh:
            values = line.strip().split('\t')
            tmps = { heads[col]: values[col]  for col in range(0,len(values))}
            data[tmps[key_title]] = tmps
    return data

def match_probe(target_data, probe_data):
    # 1. 先对数据中的坐标转数字
    for name in target_data.keys():
        target_data[name]['center'] = int(target_data[name]['center'])
        target_data[name]['start'] = int(target_data[name]['start'])
        target_data[name]['end'] = int(target_data[name]['end'])
    for probe_name in probe_data.keys():
        probe_data[probe_name]['start'] = int(probe_data[probe_name]['start'])
        probe_data[probe_name]['end'] = int(probe_data[probe_name]['end'])
        probe_data[probe_name]['5_tm'] = float(probe_data[probe_name]['5_tm'])
        probe_data[probe_name]['3_tm'] = float(probe_data[probe_name]['3_tm'])
    
    # 2. 计算距离最近的探针
    duplicate_probe = {}
    exist_probe = []  # 已有探针。针对C类探针，需要考虑上下5M范围不能有其他探针
    
    total = 0
    good = 0
    for name in target_data.keys():
        total += 1
        chrom = target_data[name]['chrom']
        center = target_data[name]['center']
        start = target_data[name]['start']
        end = target_data[name]['end']
        # 计算每一个探针的距离（染色体过滤）、Tm值偏差
        parms = {}
        for probe_name in probe_data.keys():
            # 染色体过滤
            if chrom != probe_data[probe_name]['chrom']:
                continue
            # 是否覆盖
            overlap_start = max(start, probe_data[probe_name]['start'])
            overlap_end = min(end, probe_data[probe_name]['end'])
            if overlap_end < overlap_start:
                continue
            # 探针是否已经使用
            if probe_name in duplicate_probe:
                continue
            # 不能包含snp
            if probe_data[probe_name]['SNP（MAF）'] != '.':
                continue
            # 满足条件
            parms[probe_name] = { "tm_variance": max(abs(probe_data[probe_name]['5_tm'] - bestTm), abs(probe_data[probe_name]['3_tm'] - bestTm))}
        
        
        # 取tm偏差最小值对应的探针名称
        if parms:
            good += 1
            min_tm_probe_name = sorted(parms.keys(), key=lambda x: (parms[x]['tm_variance']))[0]
            target_data[name]['distance_to_probe'] = probe_data[min_tm_probe_name]['start'] - center
            target_data[name]['probe_name'] = min_tm_probe_name
            duplicate_probe[min_tm_probe_name] = 1
            # 其他要保存的探针数据
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = probe_data[min_tm_probe_name][probe_title]
        else:
            # 赋予空值
            target_data[name]['probe_name'] = '.'
            target_data[name]['distance_to_probe'] = '.'
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = '.'
    # 人工修正指定探针序列
    good_fix = fix_by_manual(target_data, fix_file)
    
    # C类探针上下5M范围内，不能有其他探针
    skip_C_5M(target_data)
    
    
    lost = total - good - good_fix
    print(f"共需要 {total} 个探针，成功设计 {good} + {good_fix} 个探针， 已完成 {round(100 * good/total, 2)} %。 还缺少 {lost} 个探针")

def workbook_format(workbook):
    """[给workbook添加几种常用的自定义样式，方便使用]

    Args:
        workbook ([obj]): [xlsxwriter.Workbook 对象]
    """
    wb_format = {}
    # 'font_size': 12,  # 字体大小
    # 'font_color': 'black',  # 字体颜色
    # 'font_name': 'Times New Roman',  # 字体
    # 'border': 1,  # 边框大小
    # 'border_color': 'black',  # 边框颜色
    # 'bg_color': 'yellow',  # 背景色
    # 'align': 'center',  # 水平对齐
    # 'valign': 'vcenter',  # 垂直对齐

    # 标题样式
    wb_format['title_style']  = workbook.add_format({'font_size': 12, 'font_color': 'black', 'font_name': 'Times New Roman', 'border': 1, 'border_color': 'black', 'bg_color': 'yellow', 'align': 'center', 'valign': 'vcenter',})
    # 常规样式
    wb_format['normal_style'] = workbook.add_format({'font_size': 12, 'font_color': 'black', 'font_name': 'Times New Roman', 'border': 1, 'border_color': 'black', 'bg_color': 'white', 'align': 'center', 'valign': 'vcenter',})
    wb_format['normal_style_left_align'] = workbook.add_format({'font_size': 12, 'font_color': 'black', 'font_name': 'Times New Roman', 'border': 1, 'border_color': 'black', 'bg_color': 'white', 'align': 'left', 'valign': 'vcenter',})
    wb_format['normal_style_mark'] = workbook.add_format({'font_size': 12, 'font_color': 'black', 'font_name': 'Times New Roman', 'border': 1, 'border_color': 'black', 'bg_color': '#c4d79b', 'align': 'center', 'valign': 'vcenter',})
    wb_format['normal_style_red'] = workbook.add_format({'font_size': 12, 'font_color': 'black', 'font_name': 'Times New Roman', 'border': 1, 'border_color': 'black', 'bg_color': 'red', 'align': 'center', 'valign': 'vcenter',})
    return wb_format

def output(target_data, titles, output_file):
    print("分析结果：" + output_file)
    # 需要特殊颜色重点标记的列
    mark_titles = ['name', 'chrom', 'center', 'distance_to_probe', 'probe_name']
    # 缺失内容
    red_titles = ['distance_to_probe', 'probe_name']
    
    workbook = xlsxwriter.Workbook(output_file)
    wb_format = workbook_format(workbook)
    worksheet = workbook.add_worksheet('A染色体探针')
    row = 0
    for col in range(len(titles)):
        worksheet.write(row, col, titles[col], wb_format['title_style'])
    row += 1
    worksheet.freeze_panes(1, 0)  # 冻结第一行
    # worksheet.set_column(0, 40, 15) # 修改列宽
    for name in target_data.keys():
        values = [target_data[name][title] for title in titles]
        for col in range(len(titles)):
            # 颜色标记
            style = 'normal_style'
            if titles[col] in mark_titles:
                style = 'normal_style_mark'
                if titles[col] in red_titles and values[col] == '.':
                    style = 'normal_style_red'
            # 输出
            worksheet.write(row, col, values[col], wb_format[style])
        row += 1
        
    ##########
    # read me 表格
    ##########
    worksheet_readme = workbook.add_worksheet('read me')
    worksheet_readme.set_column(0, 0, 30)
    worksheet_readme.set_column(1, 1, 30)
    worksheet_readme.set_column(2, 2, 30)
    readme_info = [['sheet', 'title', 'description'],
                    ['A染色体探针', 'name', '探针名称，格式为 "chrom-[ABC][PQ][0-9]"。其中 A表示A末端，B表示着丝粒，C表示平均。P 表示染色体短臂，Q表示染色体长臂。0-9为编号。A末端类型中，数字表示距离，单位为Mb。 B着丝粒类型中，0/1 也表示距离，单位为MB。C平均类型中， 为从 10MB ~ (着丝粒位置 - 1 - 5)MB 区间内，每20M 设计一个探针，两个探针之间间隔至少5MB。'],
                    ["A染色体探针", "chrom", "染色体"],
                    ["A染色体探针", "center", "要设计探针的坐标"],
                    ["A染色体探针", "type", "参考name列的后缀"],
                    ["A染色体探针", "message", "后期人工修正时的描述信息"],
                    ["A染色体探针", "distance_to_probe", "设计的探针与center的距离（目前按照200K 范围选择最近的）"],
                    ["A染色体探针", "probe_name", "探针原始名称。从cnvplex软件提取出来的，格式： 文件名-探针名。仅供溯源参考"],
                    ["A染色体探针", "probe_chrom", "cnvplex软件设计结果：探针的染色体"],
                    ["A染色体探针", "probe_start", "cnvplex软件设计结果：探针的起始坐标"],
                    ["A染色体探针", "probe_end", "cnvplex软件设计结果：探针的终止坐标"],
                    ["A染色体探针", "probe_strand", "cnvplex软件设计结果：探针的方向"],
                    ["A染色体探针", "probe_最高同源性", "cnvplex软件注释结果：最高同源性"],
                    ["A染色体探针", "probe_CNV%", "cnvplex软件注释结果：CNV MAF"],
                    ["A染色体探针", "probe_SNP（MAF）", "cnvplex软件注释结果：SNP MAF"],
                    ["A染色体探针", "probe_5_seq", "cnvplex软件设计结果：5'探针序列"],
                    ["A染色体探针", "probe_5_length", "cnvplex软件设计结果：5'探针序列长度"],
                    ["A染色体探针", "probe_5_tm", "cnvplex软件设计结果：5'探针序列tm值"],
                    ["A染色体探针", "probe_3_seq", "cnvplex软件设计结果：3'探针序列"],
                    ["A染色体探针", "probe_3_length", "cnvplex软件设计结果：3'探针序列长度"],
                    ["A染色体探针", "probe_3_tm", "cnvplex软件设计结果：3'探针序列tm值"],
                    ["备注", "", "以上探针，是经过了以下严格条件筛选后的结果。 1：要求探针 5 3 序列连接点两侧6bp内，必须包含一个A 或 T；2：STR 限制，序列中不存在polyNm(m>5) 或polyN2/N3m(m>3) 或polyN4(m>2)或polyN5/N6/N7m(m>1)或任何8bp及以上序列重复2次及以上；3：同源性 判断，长度不超过 30；4：200bp序列的GC含量 35%~65%"]
                   ]
    row = 0
    for sheet, title, desc in readme_info:
        style = 'title_style' if row == 0 else 'normal_style_left_align'
        worksheet_readme.write(row, 0, sheet, wb_format[style])
        worksheet_readme.write(row, 1, title, wb_format[style])
        worksheet_readme.write(row, 2, desc, wb_format[style])
        row += 1
    workbook.close()


# 读入数据
target_data = read_input(target_file, 'name')
probe_data = read_input(probe_file, 'probe_name')

# 目标、探针组合
match_probe(target_data, probe_data)

# 输出
# 准备输出文件的表头：
titles_raw = []
with open(target_file, 'r') as fh:
    titles_raw = fh.readline().strip().split('\t')
titles_raw.append('distance_to_probe')
titles_raw.append('probe_name')
for probe_title in probe_titles:
    titles_raw.append('probe_' + probe_title)

# 原始内部看的数据
output(target_data, titles_raw, output_file_raw)

# 可以发给苏州的数据
# 剔除冗余数据，用于反馈给苏州
titles = [title for title in titles_raw if not title.startswith("cnvplex")]
output(target_data, titles, output_file)



# 输出选用的探针文本，方便亚区使用
with open(output_file + ".selected_probe.txt", 'w') as fh:
    titles = ['name', 'probe_chrom','probe_start', 'probe_end']
    fh.write("\t".join(titles) + "\n")
    for name in target_data.keys():
        # 有探针，且是可以保留的，那么就是设计成功的探针
        if target_data[name]['probe_chrom'] != '.':
            values = [str(target_data[name][title]) for title in titles]
            fh.write("\t".join(values) + "\n")
