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
from utils import exclude_probe, find_overlap, read_exist_probe

warnings.filterwarnings('ignore')
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout,
    level=logging.INFO
    )
log = logging.getLogger(__name__)

# 把目标信息 与 探针信息组合起来
target_file ='./noney/sub_region/sub_region_design.unsoftmask.txt'
probe_file = './probe/sub_region/final.txt'
output_file = "./probe/sub_region/sub_region_probe.xlsx"
output_file_raw = "./probe/sub_region/sub_region_probe_raw.xlsx"
exists_probes = {"core_gene": './probe/core_gene/core_gene_probe.xlsx.selected_probe.txt', "critical_region": './probe/critical_region/critical_region_probe.xlsx.selected_probe.txt', 'chromsome': './probe/chromosome/A_probe.xlsx.selected_probe.txt'}  # 核心基因、核心区域已有探针，这些区域可以不用设计了
exist_probe = read_exist_probe(exists_probes)  # 读入已有探针
bestTm = 65 # 最佳Tm 值
probe_titles = ['chrom', 'start', 'end', 'strand', '最高同源性', 'CNV%', 'SNP（MAF）', '5_seq', '5_length', '5_tm', '3_seq', '3_length', '3_tm']

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
        target_data[name]['bigzone'] = name.split('-')[0]
        target_data[name]['from'] = '亚区探针'
        target_data[name]['center'] = int(target_data[name]['center'])
        target_data[name]['start'] = int(target_data[name]['start'])
        target_data[name]['end'] = int(target_data[name]['end'])
        target_data[name]['width'] = int(target_data[name]['width'])
        target_data[name]['probe_number'] = int(target_data[name]['probe_number'])
        target_data[name]['REGION90 Sample Number'] = int(target_data[name]['REGION90 Sample Number'])
    for probe_name in probe_data.keys():
        probe_data[probe_name]['start'] = int(probe_data[probe_name]['start'])
        probe_data[probe_name]['end'] = int(probe_data[probe_name]['end'])
        probe_data[probe_name]['5_tm'] = float(probe_data[probe_name]['5_tm'])
        probe_data[probe_name]['3_tm'] = float(probe_data[probe_name]['3_tm'])
    
    # 2. 计算距离最近的探针
    duplicate_probe = {}
    
    total = 0
    good = 0
    skip = 0
    for name in sorted(target_data.keys(), key=lambda x: target_data[name]['REGION90 Sample Number'], reverse=True):  # 优先设计样本数量多的
        total += 1
        chrom = target_data[name]['chrom']
        start = target_data[name]['start']
        end = target_data[name]['end']
        center = target_data[name]['center']
        
        # 先计算要不要设计
        e_chrom, e_start, e_end = re.split('[:-]', target_data[name]['exclude_if'])
        overlaps = find_overlap({"chrom": e_chrom, "start": int(e_start), "end": int(e_end)}, exist_probe)
        if overlaps:
            skip += 1
            target_data[name]['exclude_status'] = overlaps[0]['name']
            # 赋予空值
            target_data[name]['probe_name'] = '.'
            target_data[name]['distance_to_probe'] = '.'
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = '.'
            # 当前区域不用再找探针了，直接跳过
            continue  
        else:
            target_data[name]['exclude_status'] = 'keep'
        
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
            # 满足条件
            parms[probe_name] = {"overlap_length": overlap_end - overlap_start + 1, "center": probe_data[probe_name]['start'] - center, "tm_variance": max(abs(probe_data[probe_name]['5_tm'] - bestTm), abs(probe_data[probe_name]['3_tm'] - bestTm))}
        # 取距离最小值对应的探针名称
        if parms:
            good += 1
            min_tm_probe_name = sorted(parms.keys(), key=lambda x: (parms[x]['tm_variance']))[0]
            target_data[name]['probe_name'] = min_tm_probe_name
            target_data[name]['distance_to_probe'] = parms[min_tm_probe_name]['center']
            duplicate_probe[min_tm_probe_name] = 1
            # 其他要保存的探针数据
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = probe_data[min_tm_probe_name][probe_title]
            # 保存到已有探针数据库中
            exist_probe.append({"chrom": chrom, "start": probe_data[min_tm_probe_name]['start'], "end": probe_data[min_tm_probe_name]['end'], "name": "sub_region," + name})
        else:
            # 赋予空值
            target_data[name]['probe_name'] = '.'
            target_data[name]['distance_to_probe'] = '.'
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = '.'
    lost = total - good - skip
    print(f"共需要 {total} 个探针，成功设计 {good + skip} 个探针(good={good}, skip={skip})， 已完成 {round(100 * (good + skip)/total, 2)} %。 还缺少 {lost} 个探针")

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
    mark_titles = ['name', 'gene', 'center', 'status', 'type', 'chrom', 'start', 'end', 'distance_to_probe', 'probe_name', 'exclude_status']
    # 缺失内容
    red_titles = ['distance_to_probe', 'probe_name']
    
    workbook = xlsxwriter.Workbook(output_file)
    wb_format = workbook_format(workbook)
    worksheet = workbook.add_worksheet('亚区')
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
                    ['亚区', 'name', '探针名称，格式为 "大区-chromosome_start_end-[L|R|Mn|L75|L50|L50_75|R75|R50|R50_75]"。核心区域的染色体、起始、终止。L/R 分别表示核心区域的 5% 位置和95%位置。M表示中间均匀分布的探针，M后缀还有一个数字，表示编号，从0开始。L50/L75 表示在 L-REGION50/L-REGION75 的5%的位置。R50/R75 表示在 R-REGION50/R-REGION75 的5%的位置。L50_75/R50_75 表示左右两个50 75 合并后的5%的位置。 L/R75-1 或 LR/75-2 表示75区域设置两个探针（大于2倍核心区域长度时才存在）  '],
                    ["亚区", "sub_chrom", "亚区染色体"],
                    ["亚区", "sub_start", "亚区起始"],
                    ["亚区", "sub_end", "亚区终止"],
                    ["亚区", "sub_width", "亚区宽度"],
                    ["亚区", "L-REGION50", "左侧50位置坐标"],
                    ["亚区", "L-REGION75", "左侧75位置坐标"],
                    ["亚区", "R-REGION75", "右侧50位置坐标"],
                    ["亚区", "R-REGION50", "右侧75位置坐标"],
                    ["亚区", "REGION90 Sample Number", "样本数量。后续设计探针时，优先考虑样本数量多的亚区。"],
                    ["亚区", "probe_number", "设计的探针数量"],
                    ["亚区", "interval", "探针两两之间的距离。左右两侧 50 75 区域不存在这个值"],
                    ["亚区", "center", "探针设计的核心位点。例如 5% 位点，95%位点，以及中间其他的位点"],
                    ["亚区", "chrom", "需要设计探针的区域的染色体"],
                    ["亚区", "start", "需要设计探针的区域的起始"],
                    ["亚区", "end", "需要设计探针的区域的终止"],
                    ["亚区", "width", "需要设计探针的区域的宽度"],
                    ["亚区", "exclude_if", "当该区域内存在探针时，当前name对应的探针可以取消。它是用来做后续过滤的。"],
                    ["亚区", "exclude_status", "核心基因或核心区域是否已在 exclude_if 区域内存在探针。 keep: 不存在， 其他：核心基因或亚区在exclude_if区域中的探针名"],
                    ["亚区", "distance_to_probe", "探针与center的距离"],
                    ["亚区", "probe_name", "探针原始名称。从cnvplex软件提取出来的，格式： 文件名-探针名。仅供溯源参考"],
                    ["亚区", "probe_chrom", "cnvplex软件设计结果：探针的染色体"],
                    ["亚区", "probe_start", "cnvplex软件设计结果：探针的起始坐标"],
                    ["亚区", "probe_end", "cnvplex软件设计结果：探针的终止坐标"],
                    ["亚区", "probe_strand", "cnvplex软件设计结果：探针的方向"],
                    ["亚区", "probe_最高同源性", "cnvplex软件注释结果：最高同源性"],
                    ["亚区", "probe_CNV%", "cnvplex软件注释结果：CNV MAF"],
                    ["亚区", "probe_SNP（MAF）", "cnvplex软件注释结果：SNP MAF"],
                    ["亚区", "probe_5_seq", "cnvplex软件设计结果：5'探针序列"],
                    ["亚区", "probe_5_length", "cnvplex软件设计结果：5'探针序列长度"],
                    ["亚区", "probe_5_tm", "cnvplex软件设计结果：5'探针序列tm值"],
                    ["亚区", "probe_3_seq", "cnvplex软件设计结果：3'探针序列"],
                    ["亚区", "probe_3_length", "cnvplex软件设计结果：3'探针序列长度"],
                    ["亚区", "probe_3_tm", "cnvplex软件设计结果：3'探针序列tm值"],
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
titles_raw.insert(0, 'from')  # 来源
titles_raw.insert(0, 'bigzone')  # 大区
titles_raw.append('distance_to_probe')
titles_raw.append('probe_name')
titles_raw.append('exclude_status')
for probe_title in probe_titles:
    titles_raw.append('probe_' + probe_title)
output(target_data, titles_raw, output_file_raw)

# 剔除冗余数据，用于反馈给苏州
titles = [title for title in titles_raw if not title.startswith("cnvplex")]
output(target_data, titles, output_file)