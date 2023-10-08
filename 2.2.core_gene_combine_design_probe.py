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

warnings.filterwarnings('ignore')
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout,
    level=logging.INFO
    )
log = logging.getLogger(__name__)

# 把目标信息 与 探针信息组合起来
target_file ='./noney/core_gene/core_gene_design.txt'
probe_file = './probe/core_gene/final.txt'
output_file = "./probe/core_gene/core_gene_probe.xlsx"
output_file_raw = "./probe/core_gene/core_gene_probe_raw.xlsx"

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
    
    total = 0
    good = 0
    for name in target_data.keys():
        total += 1
        
        # 可覆盖区域
        # 1. original
        need_region = []
        need_region.append({"probe_cover_source": "original", "chrom": target_data[name]['chrom'], "start": target_data[name]['start'], "end": target_data[name]['end']})
        center = int((target_data[name]['start'] + target_data[name]['end']) / 2)
        # 2. candidate exon
        if target_data[name]['candidate_exon'] != '.':
            for cnvplex in target_data[name]['candidate_exon'].split(";"):
                cnvplex_name, cnvplex_region = cnvplex.split(',')
                chrom, start, end = re.split('[:-]', cnvplex_region)
                need_region.append({"probe_cover_source": cnvplex_name.split("-")[-1], "chrom": chrom, "start": int(start), "end": int(end)})
        # 3. candidate updown
        if target_data[name]['candidate_updown'] != '.':
            for cnvplex in target_data[name]['candidate_updown'].split(";"):
                cnvplex_name, cnvplex_region = cnvplex.split(',')
                chrom, start, end = re.split('[:-]', cnvplex_region)
                need_region.append({"probe_cover_source": cnvplex_name.split("-")[-1], "chrom": chrom, "start": int(start), "end": int(end)})
        
        # 按顺序获取探针,只需要1个
        is_match = False
        for region in need_region:
            probe_cover_source, chrom, start, end = region['probe_cover_source'], region['chrom'], region['start'], region['end']
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
                parms[probe_name] = {"overlap_length": overlap_end - overlap_start + 1, "distance_to_probe": probe_data[probe_name]['start'] - center ,"tm_variance": max(abs(probe_data[probe_name]['5_tm'] - bestTm), abs(probe_data[probe_name]['3_tm'] - bestTm))}
            # 取距离最小值对应的探针名称
            if parms:
                good += 1
                min_tm_probe_name = sorted(parms.keys(), key=lambda x: (parms[x]['tm_variance']))[0]
                target_data[name]['distance_to_probe'] = parms[min_tm_probe_name]['distance_to_probe']
                target_data[name]['probe_name'] = min_tm_probe_name
                target_data[name]['probe_cover_source'] = probe_cover_source
                duplicate_probe[min_tm_probe_name] = 1
                # 其他要保存的探针数据
                for probe_title in probe_titles:
                    target_data[name]['probe_' + probe_title] = probe_data[min_tm_probe_name][probe_title]
                is_match = True
                # 只要找到了，就停止
                break
        # 是否找到？
        if not is_match:
            # 赋予空值
            target_data[name]['distance_to_probe'] = '.'
            target_data[name]['probe_name'] = '.'
            target_data[name]['probe_cover_source'] = '.'
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = '.'
    lost = total - good
    print(f"共需要 {total} 个探针，成功设计 {good} 个探针， 已完成 {round(100 * good/total, 2)} %。 还缺少 {lost} 个探针")

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
    mark_titles = ['name', 'gene', 'center', 'status', 'type', 'chrom', 'start', 'end', 'distance_to_probe', 'probe_name', 'probe_cover_source']
    # 缺失内容
    red_titles = ['distance_to_probe', 'probe_name', 'probe_cover_source']
    
    workbook = xlsxwriter.Workbook(output_file)
    wb_format = workbook_format(workbook)
    worksheet = workbook.add_worksheet('核心基因')
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
                    ['核心基因', 'name', '探针名称，格式为 "大区-基因-[123]"。1/2/3 分别表示基因的起始、中间、终止。如果转录本是+方向的，1的坐标最小。如果转录本是-方向的，1的坐标最大。'],
                    ["核心基因", "gene", "核心基因"],
                    ["核心基因", "mrna", "核心基因对应的转录本"],
                    ["核心基因", "cds_count", "转录本包含的CDS数量"],
                    ["核心基因", "strand", "转录本的转录方向"],
                    ["核心基因", "status", "核心区域与核心基因之间的覆盖情况。共6种情况。1：two_cover 起始、终止密码子都被核心区域覆盖到了。start_cover/end_cover：只有起始或者只有终止密码子被核心区域覆盖到了。inner_cover：核心基因包含了核心区域，起始、终止密码子没有被核心区域覆盖。no_critical：当前大区没有核心区域。no_cover: 核心基因与核心区域没有交集。"],
                    ["核心基因", "type", "当前探针的设计规则情况。按照设计要求来做的，共4种情况。start_codon/stop_codon： 按照起始/终止密码子覆盖的情况进行设计，设计范围为起始/终止密码子上下1k范围内。closest_CDS：按照中间点最近的外显子设计。demarcation_point：按照交界点的最近覆盖外显子设计。"],
                    ["核心基因", "middle_point", "中心点。start_codon stop_codon 下，表示起始、终止密码子坐标。closest_CDS 下，表示中间点。demarcation_point下，表示交界点。 "],
                    ["核心基因", "center", "start/end 的覆盖范围中心点。 = (start + end) / 2"],
                    ["核心基因", "chrom", "需要设计探针的区域的染色体"],
                    ["核心基因", "start", "需要设计探针的区域的起始"],
                    ["核心基因", "end", "需要设计探针的区域的终止"],
                    ["核心基因", "width", "需要设计探针的区域的宽度"],
                    ["核心基因", "candidate_exon", "备选的所有外显子，根据距离由近及远排序。当chrom:start-end 外显子区域内无法设计探针时，可在备选区域进行设计。"],
                    ["核心基因", "candidate_updown", "备选的起始、终止密码子扩展区域。当chrom:start-end 外显子区域内无法设计探针时，可在备选区域进行设计。"],
                    ["核心基因", "critical_region_raw", "当前大区的核心区域汇总"],
                    ["核心基因", "critical_region_merge", "当前大区的核心区域merge成一个区域"],
                    ["核心基因", "cds", "转录本的每一个cds区域"],
                    ["核心基因", "start_codon", "转录本的起始密码子区域"],
                    ["核心基因", "stop_codon", "转录本的终止密码子区域"],
                    ["核心基因", "distance_to_probe", "探针与  center 或candidate_center(candidate_exon的中心点) 的距离"],
                    ["核心基因", "probe_cover_source", "探针覆盖的对象：original 表示覆盖 chrom:start-end区域，candidate_CDS表示覆盖 candidate_exon 区域, expandnK，表示覆盖 candidate_updown区域"],
                    ["核心基因", "probe_name", "探针原始名称。从cnvplex软件提取出来的，格式： 文件名-探针名。仅供溯源参考"],
                    ["核心基因", "probe_chrom", "cnvplex软件设计结果：探针的染色体"],
                    ["核心基因", "probe_start", "cnvplex软件设计结果：探针的起始坐标"],
                    ["核心基因", "probe_end", "cnvplex软件设计结果：探针的终止坐标"],
                    ["核心基因", "probe_strand", "cnvplex软件设计结果：探针的方向"],
                    ["核心基因", "probe_最高同源性", "cnvplex软件注释结果：最高同源性"],
                    ["核心基因", "probe_CNV%", "cnvplex软件注释结果：CNV MAF"],
                    ["核心基因", "probe_SNP（MAF）", "cnvplex软件注释结果：SNP MAF"],
                    ["核心基因", "probe_5_seq", "cnvplex软件设计结果：5'探针序列"],
                    ["核心基因", "probe_5_length", "cnvplex软件设计结果：5'探针序列长度"],
                    ["核心基因", "probe_5_tm", "cnvplex软件设计结果：5'探针序列tm值"],
                    ["核心基因", "probe_3_seq", "cnvplex软件设计结果：3'探针序列"],
                    ["核心基因", "probe_3_length", "cnvplex软件设计结果：3'探针序列长度"],
                    ["核心基因", "probe_3_tm", "cnvplex软件设计结果：3'探针序列tm值"],
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
titles_raw.append('probe_cover_source')
titles_raw.append('probe_name')
for probe_title in probe_titles:
    titles_raw.append('probe_' + probe_title)
output(target_data, titles_raw, output_file_raw)

# 剔除冗余数据，用于反馈给苏州
titles = [title for title in titles_raw if not title.startswith("cnvplex")]
output(target_data, titles, output_file)

# 输出选用的探针文本，方便核心区、亚区使用
with open(output_file + ".selected_probe.txt", 'w') as fh:
    titles = ['name', 'probe_chrom','probe_start', 'probe_end']
    fh.write("\t".join(titles) + "\n")
    for name in target_data.keys():
        if target_data[name]['probe_chrom'] != '.':
            values = [str(target_data[name][title]) for title in titles]
            fh.write("\t".join(values) + "\n")