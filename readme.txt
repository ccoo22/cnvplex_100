1. unsoftmask.py
    从hg19基因组中，提取非小写、非N字符的区域，方便后续在这些区域设计引物
    结果文件是 target.txt

2. run_a.py
    基于 "A染色体探针" 的设计规则，找到要设计的探针区域
    输出的结果在 chromosome/final.txt 中，内容包括
    chrom	start	end	 name	       center	            type	     cnvplex_format
    染色体  起始     终止 当前探针的名字 当前探针的坐标中心点    探针的类型   方便用于cnvplex设计探针的格式

    type解释：
        1. 开头字母A/B/C,分别对应 A末端、B着丝粒、C平均，P/Q 分别表示染色体长臂、短臂
        A开头的类型中，后缀数字，表示 0、1、2、3、5MB 目标位置
        B开头的类型中，后缀数字，表示 0、1 MB 位置
        C开头的类型中，后缀数字，表示编号，第0个/第1个/...

3. probe_filter_a.py
    对cnvplex软件输出的探针excel结果解析，并按照 “探针分布设计规则” 中的要求，对探针进行过滤，输出合格的探针汇总结果
    输入在 probe_a 
          所有的excel表格，可以直接放在该目录下，软件会全部读入
    输出在 probe_a_filter/final.txt

4. run_a_unsoftmask.py
   对 "A染色体探针" 的目标区域处理，提取探针区域内