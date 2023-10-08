#!/home/genesky/software/python/3.9.4/bin/python3
import logging
import multiprocessing
import os
import re
import sys

from utils import design_cnvplex, fix_design_center

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout,
    level=logging.INFO
    )
log = logging.getLogger(__name__)



def create_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

def read_genome_length(fai):
    genome_length = {}
    with open(fai, 'r') as fh:
        for line in fh:
            values = line.strip().split()
            genome_length[values[0]] = int(values[1])
    return genome_length

def read_centromere_bed(bed):
    bed_region = {}
    with open(bed, 'r') as fh:
        for line in fh:
            chrom, start, end = line.strip().split()
            bed_region[chrom] = {'start': int(start), 'end': int(end)}
    return bed_region

# 要分析的染色体
chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
# P末端不需要设计的染色体
chroms_skip = ["chr13", "chr14", "chr15","chr21", "chr22"]
# 基因组文件
genome_fa = '/home/genesky/database_new/ucsc/fasta/hg19/samtools_index/hg19.fa'
# 基因组长度文件
genome_fai = '/home/genesky/database_new/ucsc/fasta/hg19/samtools_index/hg19.fa.fai'
genome_length = read_genome_length(genome_fai)
# 质粒bed文件
centromere_bed = "./centromere.bed"
centromere_region = read_centromere_bed(centromere_bed)
# ATCG 覆盖区域
target_bed = "./target.bed"
# samtools 软件
samtools = '/home/genesky/software/samtools/1.17/bin/samtools'
# 结果输出目录
output_dir = "./chromosome"

# 后期人工修正
fix_file = "./buchong/chromsome.txt"

create_dir(output_dir)
steps = [10, 50, 100, 200] # 每次跨越步长

def my_func(chrom):
    # log.info(f"[process] chromosome {chrom} ")
    chrom_length = genome_length[chrom]
    centromere_start, centromere_end = centromere_region[chrom]['start'], centromere_region[chrom]['end']
    # step = 100*1000  
    KB1 = 1000  # 1k大小
    MB1 = KB1 * KB1 # 1mb 大小
    
    interval_min = 5 * MB1 # 染色体臂上的两个探针之间最小间隔 5M
    interval_max = 20 * MB1 # 染色体臂上的两个探针之间最大间隔 20M
    # 要设计的区域汇总
    designed_region = {}

    ### 1. 设计A末端
    moduan = [0, 1, 2, 3, 5]
    if chrom not in chroms_skip:  # 指定的染色体不做
        # 1.1 P末端
        for p in moduan:
            type = "AP" + str(p)
            name = chrom + "-" + type
            center = p * MB1
            designed_region[name] = {'name': name}
            designed_region[name]['chrom'] = chrom
            designed_region[name]['center'] = center if(p != 0) else 1
            designed_region[name]['type'] = type
            designed_region[name]['expand'] = 200 * KB1  # 允许上下扩展的范围
            
            # # 生成每种expand下的坐标
            # for step in steps:
            #     start = center - step * KB1
            #     end = center + step * KB1
            #     if start < 0:
            #         start = 1
            #     designed_region[name]['expand_' + str(step) + "K"] = name + "," + chrom + ":" + str(start) + "-" + str(end)
    # 1.2 q末端
    for q in moduan:
        type = "AQ" + str(q)
        name = chrom + "-" + type
        center = chrom_length - q * MB1
        designed_region[name] = {'name': name}
        designed_region[name]['chrom'] = chrom
        designed_region[name]['center'] = center
        designed_region[name]['type'] = type
        designed_region[name]['expand'] = 200 * KB1  # 允许上下扩展的范围
        # 生成每种expand下的坐标
        # for step in steps:
        #     start = center - step * KB1
        #     end = center + step * KB1
        #     if end > chrom_length:
        #         end = chrom_length
        #     designed_region[name]['expand_' + str(step) + "K"] = name + "," + chrom + ":" + str(start) + "-" + str(end)


    ### 2. 着丝粒
    centers = [centromere_start - MB1, centromere_start, centromere_end, centromere_end + MB1]
    types = ["BP1", "BP0", "BQ0", "BQ1"]
    for type, center in zip(types, centers):
        name = chrom + "-" + type
        # 特殊染色体的P端不用设计
        if chrom in chroms_skip and name.startswith('BP'):
            continue
        designed_region[name] = {'name': name}
        designed_region[name]['chrom'] = chrom
        designed_region[name]['center'] = center
        designed_region[name]['type'] = type
        designed_region[name]['expand'] = 200 * KB1  # 允许上下扩展的范围
        # 生成每种expand下的坐标
        # for step in steps:
        #     start = center - step * KB1
        #     end = center + step * KB1
        #     designed_region[name]['expand_' + str(step) + "K"] = name + "," + chrom + ":" + str(start) + "-" + str(end)
            
    ### 3. C平均
    # 3.1 P端
    if chrom not in chroms_skip:
        begin = 10 * MB1  # 起点为10MB
        end = centromere_start - 6 *MB1  # 终点
        poses = list(range(begin, end, interval_max))
        poses_select = []
        # print(chrom, begin, end, "P", end - begin, poses)
        # 不足
        if len(poses) == 0:
            pass
        elif len(poses) == 1:
            # 只有1个，需要判断
            if (end - begin) > interval_min:
                poses_select = [begin, end]
            else:
                poses_select = [(begin + end) / 2]
        else:
            poses_select = poses
        
        for index, center in enumerate(poses_select):
            type = "CP" + str(index + 1)
            name = chrom + "-" + type
            designed_region[name] = {'name': name}
            designed_region[name]['chrom'] = chrom
            designed_region[name]['center'] = center
            designed_region[name]['type'] = type
            designed_region[name]['expand'] = MB1  # 允许上下扩展的范围
            # # 生成每种expand下的坐标
            # for step in steps:
            #     start = center - step * KB1
            #     end = center + step * KB1
            #     designed_region[name]['expand_' + str(step) + "K"] = name + "," + chrom + ":" + str(start) + "-" + str(end)
    
    # 3.2 Q端
    begin = centromere_end + 6 * MB1  # 起点
    end = chrom_length - 10 * MB1  # 终点
    poses = list(range(begin, end, interval_max))
    # print(chrom, begin, end, "P", end-begin, poses)
    # 不足
    if len(poses) == 0:
        pass
    elif len(poses) == 1:
        # 只有1个，需要判断
        if (end - begin) > interval_min:
            poses_select = [begin, end]
        else:
            poses_select = [(begin + end) / 2]
    else:
        poses_select = poses
    
    for index, center in enumerate(poses_select):
        type = "CQ" + str(index + 1)
        name = chrom + "-" + type
        designed_region[name] = {'name': name}
        designed_region[name]['chrom'] = chrom
        designed_region[name]['center'] = center
        designed_region[name]['type'] = type
        designed_region[name]['expand'] = MB1  # 允许上下扩展的范围
        # # 生成每种expand下的坐标
        # for step in steps:
        #     start = center - step * KB1
        #     end = center + step * KB1
        #     designed_region[name]['expand_' + str(step) + "K"] = name + "," + chrom + ":" + str(start) + "-" + str(end)
    
    
    # 补充修正
    fix_design_center(designed_region, fix_file )
    
    # 允许设计探针的区域
    for name in designed_region.keys():
        start = designed_region[name]['center'] - designed_region[name]['expand']
        end = designed_region[name]['center'] + designed_region[name]['expand']
        if start <= 0:
            start = 1
        if end > chrom_length:
            end = chrom_length
        designed_region[name]['start'] = start
        designed_region[name]['end'] = end
        designed_region[name]['width'] = end - start + 1
        if 'message' not in designed_region[name]:
            designed_region[name]['message'] = '.'

    
    # 制作cnvplex 格式数据
    design_cnvplex(designed_region)
    
    ### 
    # 区域信息输出
    chrome_res = os.path.join(output_dir, chrom + ".txt")
    titles = ['name', 'center', 'type', 'chrom', 'start', 'end', 'width', 'message', 'cnvplex_1K', 'cnvplex_10K', 'cnvplex_100K', 'cnvplex_200K']

    with open(chrome_res, 'w') as fh:
        for name in designed_region.keys():
            values = [str(designed_region[name][title])  for title in titles]
            fh.write("\t".join(values) + "\n")

pool = multiprocessing.Pool(processes=10)
for chrom in chroms:
    pool.apply_async(func=my_func, args=(chrom, ), error_callback = lambda e:print(e.__cause__))
pool.close()
pool.join()


# 汇总
titles = ['name', 'center', 'type', 'chrom', 'start', 'end', 'width', 'message', 'cnvplex_1K', 'cnvplex_10K', 'cnvplex_100K', 'cnvplex_200K']

with open(f"{output_dir}/final.txt", "w") as fh:
    fh.write("\t".join(titles) + "\n")

final_file = f"{output_dir}/final.txt"
for chrom in chroms:
    os.system(f"cat {output_dir}/{chrom}.txt >> {final_file}")

print("结果文件：" + final_file)