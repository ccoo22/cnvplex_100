#!/home/genesky/software/python/3.8.12/bin/python3
import os
import re


# 解析gtf注释列，拆分成字典
def parse_gtf_anno(annotation_string):
    annotation_dict = {}
    for anno in annotation_string.split('";'):
        if not re.search('\w', anno):
            continue
        name, value = anno.strip().split(' "', 1)
        annotation_dict[name] = value
    return annotation_dict


# 解析gtf文件，提取目标转录本
def parse_gtf(gtf):
    print(f"解析GTF {gtf}")
    features = ['exon', 'CDS', 'start_codon', 'stop_codon']
    gtf_info = {}
    with open(gtf, 'r') as fh:
        for line in fh:
            if re.match('#', line):
                continue
            # 提取gtf一行中的信息
            chrom, source, feature, start, end, score, strand, phread, annotation_string = line.strip().split('\t')
            if chrom.startswith('NW') or chrom.startswith('NT'):
                continue
            start, end = int(start), int(end)
            annotation_dict = parse_gtf_anno(annotation_string)
            gene_id = annotation_dict['gene_id']
            transcript_id = annotation_dict['transcript_id']
            if transcript_id == '':
                continue
            mrna = transcript_id.split('.')[0]
            if mrna not in gtf_info:
                gtf_info[mrna] = {}
                
            # 读入转录本的exon/start_codon/stop_codon/CDS 信息
            if feature in features:
                if feature not in gtf_info[mrna]:
                    gtf_info[mrna][feature] = []
                feature_number = annotation_dict.get('exon_number', '.')
                gtf_info[mrna][feature].append({'chrom': chrom, 'start': start, 'end': end, 'strand': strand, 'feature': feature, 'feature_number': feature_number})
    return gtf_info
    

def read_core_gene(file):
    '''读取核心基因文件'''
    core_gene = {}
    with open(file, 'r') as fh:
        fh.readline()
        for line in fh:
            gene, trans, oldname = line.strip().split('\t')
            if oldname != '.':
                gene = oldname
            trans = trans.split('.')[0]
            core_gene[gene] = trans
            
    return core_gene


def read_core_gene_trans(file):
    '''读取核心基因的转录本区域信息'''
    core_gene_region = {}
    with open(file, 'r') as fh:
        heads = fh.readline().strip().split('\t')
        for line in fh:
            values = line.strip().split('\t')
            tmps = { heads[col]: values[col]  for col in range(0,len(values))}
            feature = tmps['feature']
            if tmps['gene'] not in core_gene_region:
                core_gene_region[tmps['gene']] = {'mrna': tmps['trans']}
            if feature not in core_gene_region[tmps['gene']]:
                core_gene_region[tmps['gene']][feature] = []
            core_gene_region[tmps['gene']][feature].append({'chrom': tmps['chrom'], 'start': int(tmps['start']), 'end': int(tmps['end']), 'strand': tmps['strand'], 'feature': tmps['feature'], 'feature_number': tmps['feature_number'], "name": tmps['feature'] + tmps['feature_number']})
    return core_gene_region

def find_overlap(one_obj, others_obj):
    '''计算一个区域与其他区域是否有交集，以及交集的坐标'''
    # chrom/start/end
    overlaps = []
    for index, other_obj in enumerate(others_obj):
        if one_obj['chrom'] == other_obj['chrom']:
            start = max(one_obj['start'], other_obj['start'])
            end = min(one_obj['end'], other_obj['end'])
            name = '.'
            if 'name' in other_obj:
                name = other_obj['name']
            # 存在交集
            if end >= start:
                overlaps.append({"chrom": one_obj['chrom'], "start": start, "end": end, "overlap_length": end - start + 1, 'index': index, 'name': name})  # index ： 与哪个other_obj 存在交集
    return overlaps

def link_region(objs):
    '''把多个独立的区域，连接成一个大区域: 染色体需要一样'''
    chrom = ''
    start = 0
    end = 0
    for index, obj in enumerate(objs):
        if index == 0:
            chrom = obj['chrom']
            start = obj['start']
            end = obj['end']
        else:
            if chrom != obj['chrom']:
                print("[Error] link_region 中的所有objs 的染色体必须一致，否则，仅合并与第一个染色体一致的区域")
                print(objs)
                continue
            if start > obj['start']:
                start = obj['start']
            if end < obj['end']:
                end = obj['end']
    return {"chrom": chrom, "start": start, "end": end}

def find_closest_region(one_obj, others_obj, need = 1):
    '''找到与目标区域中心点最近的一个区域'''
    center = int((one_obj['start'] + one_obj['end']) / 2)
    distance = {}
    for index, other_obj in enumerate(others_obj):
        if one_obj['chrom'] != other_obj['chrom']:
            continue
        # print(center, min_distance, min_index, index, other_obj,)
        distance[index] = min(abs(center - other_obj['start']), abs(center - other_obj['end']))
    ordered_index = sorted(distance.keys(), key=lambda index: distance[index])
    if need == 1:
        return others_obj[ordered_index[0]]
    if need == 2:
        one = others_obj[ordered_index[0]]
        two = []
        if len(others_obj) > 1:
            for index in ordered_index[1:]:
                tmp = others_obj[ordered_index[index]]
                two.append(tmp['name'] + "," + tmp['chrom'] + ":" + str(tmp['start']) + "-" + str(tmp['end']))
        if two:
            two = ";".join(two)
        else:
            two = '.'
        return one, two


def read_noney(file):
    '''读取 非Y染色体微缺失微重复区域 信息'''
    noney = {}
    number_titles = ["L-REGION50", "L-REGION75", "REGION90-P1", "REGION90-P2", "R-REGION75", "R-REGION50", "REGION90 Sample Number"]
    row = 0
    with open(file, 'r') as fh:
        heads = fh.readline().strip().split('\t')
        for line in fh:
            row += 1
            values = line.strip().split('\t')
            tmps = {}
            for col, head in enumerate(heads):
                value = values[col] if (col < len(values)) else ''
                tmps[head] = value
                # tmps = { heads[col]: values[col]  for col in range(0,len(values))}
            # 数据做修正
            # 1. 染色体
            tmps['Chr'] = 'chr' + tmps['Chr']
            
            # 2. 核心区域坐标
            critical = tmps['Critical Region(s)']
            if critical != '': 
                critical_start, critical_end = critical.split('-')
                critical_start = int(critical_start.replace(',',''))
                critical_end = int(critical_end.replace(',',''))
                tmps['critical_start'] = critical_start
                tmps['critical_end'] = critical_end
            else:
                tmps['critical_start'] = 0
                tmps['critical_end'] = 0
            
            # 3. 核心基因
            if tmps['核心区关键基因'] == "":
                tmps["核心区关键基因"] = []
            else:
                tmps['核心区关键基因'] = tmps['核心区关键基因'].replace(" ", '').split(',')
            
            # 4. 亚区坐标
            for number_title in number_titles:
                if not re.search('\d', tmps[number_title]):
                    tmps[number_title] = 0
                else:
                    if re.search('\.', tmps[number_title]):
                        # 浮点型数据
                        tmps[number_title] = int(float(tmps[number_title].replace(',','')))
                    else:
                        tmps[number_title] = int(tmps[number_title].replace(',',''))
            
            # 数据保存
            bigzone = tmps["大区"]
            if bigzone not in noney:
                noney[bigzone] = {}
            noney[bigzone]['row' + str(row)] = tmps
    return noney


def region_probe_design(chrom, start, end, expand_rate = 0.1):
    '''根据大区核心区域探针设计的规则，设计探针'''
    MB1 = 1000 * 1000 # 1Mb
    length = end - start + 1
    # 探针数量
    probe_number = 0
    # 两侧位点
    besides = []
    # 中间均匀分布位点
    middles = []
    # 探针之间的间距
    interval = 0
    if length <= 500:
        probe_number = 2
        besides.append( start + int(length * 0.05))
        besides.append( start + int(length * 0.95))
        interval = besides[1] - besides[0]
    elif length <= MB1:
        probe_number = 3
        besides.append(start + int(length * 0.05))
        middles.append(start + int(length * 0.5))
        besides.append(start + int(length * 0.95))
        interval = int((besides[1] - besides[0])/ 2)
    else:
        need_probe =  len(range(int(0.05 * length), int(length * 0.95), MB1)) + 1  # 开头 5%， 截尾 5% 固定有探针，中间至少1/MB
        if need_probe <= 4:
            # 至少4个探针
            probe_number = 4
        elif need_probe <= 6:
            probe_number = need_probe
        elif need_probe > 6:
            # 最多6个探针
            probe_number = 6
        besides.append( start + int(length * 0.05))
        besides.append( start + int(length * 0.95))
        # 中间均匀分布剩下的探针
        interval = int((besides[1] - besides[0]) / (probe_number - 1))
        for i in range(probe_number - 2):
            middles.append( start + int(length * 0.05) + interval * (i +1))
    # 平均探针间隔距离
    interval2 = length / (probe_number - 1)
    # 创建要设计探针的区域
    tmps = [
        {"surfix": "L", "interval": interval,"chrom": chrom, "center": besides[0], "expand": int(length * 0.05)},
        {"surfix": "R", "interval": interval,"chrom": chrom, "center": besides[1], "expand": int(length * 0.05)}
    ]
    for index, position in enumerate(middles):
        tmps.append({"surfix": "M" + str(index), "interval": interval,"chrom": chrom, "center": position, "expand": int(interval * expand_rate)})
    # start/end 暂定为 interval 的 10%
    target_regions = []
    for tmp in tmps:
        # expect_expand = max_expand if max_expand < 100 * 1000 else 100 * 1000  # 理想扩展
        tmp['start'] = tmp['center'] - tmp['expand']
        tmp['end'] = tmp['center'] + tmp['expand']
        tmp['width'] = tmp['end'] - tmp['start'] + 1
        tmp['probe_number'] = probe_number
        # 如果该exclude_if 区域内已经有探针了，当前center就不用设计
        exclude_start = tmp['center'] - int(0.45 * interval2)
        if exclude_start < 1:
            exclude_start = 1
        exclude_end = tmp['center'] + int(0.45 * interval2)
        tmp['exclude_if'] = chrom + ":" + str(exclude_start) + "-" + str(exclude_end)
        target_regions.append(tmp)
    return target_regions


def create_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
        


def design_75_1(name, designed_region, side, info75):
    '''为75%区域设计1个探针'''
    width_005 = int(info75['width'] * 0.05)  # 5% 区域长度
    width_05 = int(info75['width'] * 0.5)  # 50% 区域长度
    designed_region[name] = {"side": side}
    if side == "L":
        designed_region[name]['center'] = info75['start'] + width_005 
        designed_region[name]['start'] = designed_region[name]['center']
        designed_region[name]['end'] = info75['start'] + width_05
        # 如果某个区域已经有探针了，这个探针可以不用设计了
        designed_region[name]['exclude_if'] = info75["chrom"] + ":" + str(info75['start']) + "-" + str(designed_region[name]['end'])
    else:
        designed_region[name]['center'] = info75['start'] - width_005 
        designed_region[name]['start'] = info75['start'] - width_05
        designed_region[name]['end'] = designed_region[name]['center']
        designed_region[name]['exclude_if'] = info75["chrom"] + ":" + str(designed_region[name]['start']) + "-" + str(info75['start'])
    designed_region[name]['width'] = designed_region[name]['end'] - designed_region[name]['start'] + 1
    designed_region[name]['probe_number'] = 1

def design_75_2(prefix, designed_region, side, info75):
    '''为75%区域设计2个探针'''
    width_005 = int(info75['width'] * 0.05)  # 5% 区域长度
    width_033 = int(info75['width'] * 0.3333333)  # 33.333% 区域长度
    
    if side == "L":
        # 1. 5% 位置设计一个
        name1 = prefix + "-1"
        designed_region[name1] = {"side": side}
        designed_region[name1 ]['center'] = info75['start'] + width_005 
        designed_region[name1]['start'] = designed_region[name1]['center']
        designed_region[name1]['end'] = info75['start'] + width_033
        # 如果某个区域已经有探针了，这个探针可以不用设计了. special_2_31: 如果区域内包含2个以上探针，且探针之间的距离大于区域的三分之一，则该探针可以不用设计
        designed_region[name1]['exclude_if'] = info75["chrom"] + ":" + str(info75['start']) + "-" + str(designed_region[name1]['end']) + ",special_2_31," + info75['chrom'] + ":" + str(info75['start']) + "-" + str(info75['start'] + info75['width'])
        designed_region[name1]['width'] = designed_region[name1]['end'] - designed_region[name1]['start'] + 1
        designed_region[name1]['probe_number'] = 2
        
        # 2. 66.66 % 位置设计一个
        name2 = prefix + "-2"
        designed_region[name2] = {"side": side}
        designed_region[name2]['center'] = info75['start'] + width_033 * 2
        designed_region[name2]['start'] = designed_region[name2]['center']
        designed_region[name2]['end'] = info75['start'] + info75['width']
        designed_region[name2]['exclude_if'] = info75["chrom"] + ":" + str(info75['start'] + width_033) + "-" + str(designed_region[name2]['end']) + ",special_2_31," + info75['chrom'] + ":" + str(info75['start']) + "-" + str(info75['start'] + info75['width'])
        designed_region[name2]['width'] = designed_region[name2]['end'] - designed_region[name2]['start'] + 1
        designed_region[name2]['probe_number'] = 2
        
    else:
        # 1. 5% 位置设计一个
        name1 = prefix + "-1"
        designed_region[name1] = {"side": side}
        designed_region[name1 ]['center'] = info75['start'] - width_005 
        designed_region[name1]['start'] = info75['start'] - width_033
        designed_region[name1]['end'] = designed_region[name1]['center']
        # 如果某个区域已经有探针了，这个探针可以不用设计了. special_2_31: 如果区域内包含2个以上探针，且探针之间的距离大于区域的三分之一，则该探针可以不用设计
        designed_region[name1]['exclude_if'] = info75["chrom"] + ":" + str(designed_region[name1]['start']) + "-" + str(info75['start']) + ",special_2_31," + info75['chrom'] + ":" + str(info75['start'] - info75['width']) + "-" + str(info75['start'])
        designed_region[name1]['width'] = designed_region[name1]['end'] - designed_region[name1]['start'] + 1
        designed_region[name1]['probe_number'] = 2
        
        # 2. 66.66 % 位置设计一个
        name2 = prefix + "-2"
        designed_region[name2] = {"side": side}
        designed_region[name2]['center'] = info75['start'] - width_033 * 2
        designed_region[name2]['start'] =  info75['start'] - info75['width']
        designed_region[name2]['end'] = designed_region[name2]['center']
        designed_region[name2]['exclude_if'] = info75["chrom"] + ":" + str(designed_region[name2]['start']) + "-" + str(info75['start'] - width_033) + ",special_2_31," + info75['chrom'] + ":" + str(info75['start'] - info75['width']) + "-" + str(info75['start'])
        designed_region[name2]['width'] = designed_region[name2]['end'] - designed_region[name2]['start'] + 1
        designed_region[name2]['probe_number'] = 2


def design_cnvplex(designed_region):
    K1 = 1000
    K10 = 10 * K1
    K50 = 50 * K1
    K100 = 100 * K1
    K200 = 200 * K1
    # 1K 10K  100K 200K 三个单位
    for name in designed_region.keys():
        chrom = designed_region[name]['chrom']
        if not re.search("[LR]\d", name):
            # 亚区
            # 0. 1K范围内
            if (designed_region[name]['width'] / 2) <= K1:
                designed_region[name]['cnvplex_1K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
            else:
                new_start = designed_region[name]['center'] - K1
                if new_start < 1:
                    new_start = 1
                designed_region[name]['cnvplex_1K'] = name + "," + chrom + ":" + str(new_start) + "-" + str(designed_region[name]['center'] + K1)
            # 1. 10K范围内
            if (designed_region[name]['width'] / 2) <= K1:
                designed_region[name]['cnvplex_10K'] = '.'
            elif (designed_region[name]['width'] / 2) <= K10:
                designed_region[name]['cnvplex_10K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
            else:
                new_start = designed_region[name]['center'] - K10
                if new_start < 1:
                    new_start = 1
                designed_region[name]['cnvplex_10K'] = name + "," + chrom + ":" + str(new_start) + "-" + str(designed_region[name]['center'] + K10)
            # 1. 50K范围内
            if (designed_region[name]['width'] / 2) <= K10:
                designed_region[name]['cnvplex_50K'] = '.'
            elif (designed_region[name]['width'] / 2) <= K50:
                designed_region[name]['cnvplex_50K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
            else:
                new_start = designed_region[name]['center'] - K50
                if new_start < 1:
                    new_start = 1
                designed_region[name]['cnvplex_50K'] = name + "," + chrom + ":" + str(new_start) + "-" + str(designed_region[name]['center'] + K50)
            # 2. 100k范围
            if (designed_region[name]['width'] / 2) <= K50:
                designed_region[name]['cnvplex_100K'] = '.'
            elif (designed_region[name]['width'] / 2) <= K100:
                designed_region[name]['cnvplex_100K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
            else:
                new_start = designed_region[name]['center'] - K100
                if new_start < 1:
                    new_start = 1
                designed_region[name]['cnvplex_100K'] = name + "," + chrom + ":" + str(new_start) + "-" + str(designed_region[name]['center'] + K100)
            # 3. 200K范围
            if (designed_region[name]['width'] / 2) <= K100:
                designed_region[name]['cnvplex_200K'] = '.'
            elif (designed_region[name]['width'] / 2) <= K200:
                designed_region[name]['cnvplex_200K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
            else:
                new_start = designed_region[name]['center'] - K200
                if new_start < 1:
                    new_start = 1
                designed_region[name]['cnvplex_200K'] = name + "," + chrom + ":" + str(new_start) + "-" + str(designed_region[name]['center'] + K200)
        else:
            # 亚区左侧
            if re.search("L\d", name):
                # 0. 1K范围内
                if (designed_region[name]['width']) <= K1:
                    designed_region[name]['cnvplex_1K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_1K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['start'] + K1)
                # 1. 10K范围内
                if (designed_region[name]['width']) <= K1:
                    designed_region[name]['cnvplex_10K'] = '.'
                elif (designed_region[name]['width']) <= K10:
                    designed_region[name]['cnvplex_10K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_10K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['start'] + K10)
                # 1. 50K范围内
                if (designed_region[name]['width']) <= K10:
                    designed_region[name]['cnvplex_50K'] = '.'
                elif (designed_region[name]['width']) <= K50:
                    designed_region[name]['cnvplex_50K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_50K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['start'] + K50)
                # 2. 100k范围
                if (designed_region[name]['width']) <= K50:
                    designed_region[name]['cnvplex_100K'] = '.'
                elif (designed_region[name]['width']) <= K100:
                    designed_region[name]['cnvplex_100K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_100K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['start'] + K100)
                # 3. 200K范围
                if (designed_region[name]['width']) <= K100:
                    designed_region[name]['cnvplex_200K'] = '.'
                elif (designed_region[name]['width']) <= K200:
                    designed_region[name]['cnvplex_200K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_200K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['start'] + K200)
            else:
                # 亚区右侧
                # 0. 1K范围内
                if (designed_region[name]['width']) <= K1:
                    designed_region[name]['cnvplex_1K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_1K'] = name + "," + chrom + ":" + str(designed_region[name]['end'] - K1) + "-" + str(designed_region[name]['end'])
                # 1. 10K范围内
                if (designed_region[name]['width']) <= K1:
                    designed_region[name]['cnvplex_10K'] = '.'
                elif (designed_region[name]['width']) <= K10:
                    designed_region[name]['cnvplex_10K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_10K'] = name + "," + chrom + ":" + str(designed_region[name]['end'] - K10) + "-" + str(designed_region[name]['end'])
                # 1. 50K范围内
                if (designed_region[name]['width']) <= K10:
                    designed_region[name]['cnvplex_50K'] = '.'
                elif (designed_region[name]['width']) <= K50:
                    designed_region[name]['cnvplex_50K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_50K'] = name + "," + chrom + ":" + str(designed_region[name]['end'] - K50) + "-" + str(designed_region[name]['end'])
                # 2. 100k范围
                if (designed_region[name]['width']) <= K10:
                    designed_region[name]['cnvplex_100K'] = '.'
                elif (designed_region[name]['width']) <= K100:
                    designed_region[name]['cnvplex_100K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_100K'] = name + "," + chrom + ":" + str(designed_region[name]['end'] - K100) + "-" + str(designed_region[name]['end'])
                # 3. 200K范围
                if (designed_region[name]['width']) <= K100:
                    designed_region[name]['cnvplex_200K'] = '.'
                elif (designed_region[name]['width']) <= K200:
                    designed_region[name]['cnvplex_200K'] = name + "," + chrom + ":" + str(designed_region[name]['start']) + "-" + str(designed_region[name]['end'])
                else:
                    designed_region[name]['cnvplex_200K'] = name + "," + chrom + ":" + str(designed_region[name]['end'] - K200) + "-" + str(designed_region[name]['end'])


def read_exist_probe(objs):
    # 1. 读入已有探针
    exist_probe = []
    for source_name in objs.keys():
        with open(objs[source_name], 'r') as fh:
            fh.readline()
            for line in fh:
                name, chrom, start, end = line.strip().split('\t')
                start = int(start)
                end = int(start)
                exist_probe.append({"chrom": chrom, "start": start, "end": end, "name": source_name + "," + name})
    return exist_probe

def exclude_probe(target_data, exists_probes):
    # 1. 读入已有探针
    exist_probe = read_exist_probe(exists_probes)

    # 2. 根据 exclude_if 的条件，判断是否探针已有
    for name in target_data.keys():
        e_chrom, e_start, e_end = re.split('[:-]', target_data[name]['exclude_if'])
        overlaps = find_overlap({"chrom": e_chrom, "start": int(e_start), "end": int(e_end)}, exist_probe)
        # 有重叠：给出剔除原因
        if overlaps:
            target_data[name]['exclude_status'] = overlaps[0]['name']
        # 无重叠：保留
        else:
            target_data[name]['exclude_status'] = 'keep'
            

def read_fix_data(file):
    fix_data = {}
    with open(file, 'r') as fh:
        fh.readline()
        for line in fh:
            name, message, center, seq5, seq3, probe_region = line.strip().split('\t')
            if center != ".":
                center = int(float(center.replace('MB', '')) * 1000 * 1000)
            if probe_region != '.':
                chrom, start, end = re.split('[:-]', probe_region)
                probe_region = {"chrom": chrom, 'start': int(start), 'end': int(end), 'strand': '.', '最高同源性': '.', 'CNV%': '.', 'SNP（MAF）': '.', '5_seq': seq5, '5_length': len(seq5), '5_tm': '.', '3_seq': seq3, '3_length': len(seq3), '3_tm': '.' }
            fix_data[name] = {
                "message": message,
                "center": center,
                'probe_region': probe_region
            }
    return fix_data

# 重新定义center
def fix_design_center(designed_region, fix_file):
    fix_data = read_fix_data(fix_file)
    for name in fix_data.keys():
        if name in designed_region:
            designed_region[name]['message'] = fix_data[name]['message']
            # 如果有重新定义center，则用新的center
            if fix_data[name]['center'] != '.':
                designed_region[name]['center'] = fix_data[name]['center']

# 人工修正探针
def fix_by_manual(target_data, fix_file):
    fix_data = read_fix_data(fix_file)
    probe_titles = ['chrom', 'start', 'end', 'strand', '最高同源性', 'CNV%', 'SNP（MAF）', '5_seq', '5_length', '5_tm', '3_seq', '3_length', '3_tm']
    good = 0  # 修复的数量
    
    for name in fix_data.keys():
        if fix_data[name]['probe_region'] != '.':
            if target_data[name]['probe_name'] == '.':
                good += 1
            target_data[name]['probe_name'] = '人工补充'
            target_data[name]['distance_to_probe'] = fix_data[name]['probe_region']['start'] - target_data[name]['center']
            for probe_title in probe_titles:
                target_data[name]['probe_' + probe_title] = fix_data[name]['probe_region'][probe_title]
    return good

def skip_C_5M(target_data):
    AB_probes = []
    for name in target_data.keys():
        if target_data[name]['type'].startswith('A') or target_data[name]['type'].startswith('B'):
            if target_data[name]['probe_chrom'] != '.':
                AB_probes.append({"chrom": target_data[name]['probe_chrom'], "start": target_data[name]['probe_start'], "end": target_data[name]['probe_end'], "name": name})
    # 对C类探针过滤
    M5 = 5 * 1000 * 1000 - 400 * 1000
    match = 0
    for name in target_data.keys():
        if target_data[name]['type'].startswith('C') and target_data[name]['probe_chrom'] != '.':
            overlaps = find_overlap({"chrom": target_data[name]['probe_chrom'], "start": int(target_data[name]['probe_start'] - M5), "end": int(target_data[name]['probe_end'] + M5)}, AB_probes)
            if overlaps:
                target_data[name]['message'] += ' skip because of ' + overlaps[0]['name']
                match + 1