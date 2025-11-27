import os
from pathlib import Path
import pandas as pd

# --------------------------
# 配置加载与目录定义（绝对路径+自动创建）
# --------------------------
configfile: "config.yaml"

# 基础目录配置（从config读取，支持自定义）
TOOL_DIR = Path(config["tool_dir"]).resolve()  # 工具脚本目录（含get_anno/、pipelines/）
INPUT_DIR = Path(config["input_dir"]).resolve()  # 原始输入Fastq目录
OUTPUT_ROOT = Path(config["output_root"]).resolve()  # 所有输出根目录
Ref = Path(config["reference_dir"]).resolve()

# 细分输出目录（自动创建，避免手动mkdir）
TEST_DIR = OUTPUT_ROOT / "test_dir"  # 修剪/去重后的Fastq目录
ANNO_DIR = OUTPUT_ROOT / "anno_files"
GLORI_OUTPUT = OUTPUT_ROOT / "new_annotated_Acites_output"  # GLORI最终输出
LOG_DIR = OUTPUT_ROOT / "logs"  # 所有规则日志目录

# 自动创建输出目录（ANNO_DIR需手动确认文件存在，此处仅创建其他目录）
for dir_path in [TEST_DIR, ANNO_DIR, GLORI_OUTPUT, LOG_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# 验证提前下载的关键文件是否存在（提前报错，避免后续规则失败）
required_files = [
    Ref / "GCF_000001405.39_GRCh38.p13_genomic.gtf",  # 提前下载的GTF（未压缩）
    Ref / "GCF_000001405.39_GRCh38.p13_assembly_report.txt",  # 提前下载的装配报告
    Ref / "GCF_000001405.39_GRCh38.p13_rna.fna",       # 提前下载的RNA序列
    Ref / "hg38.fa"                                   # 提前下载的基因组序列
]
for file in required_files:
    if not file.exists():
        raise FileNotFoundError(f"提前下载的文件不存在：{file}\n请确认文件放在 {ANNO_DIR} 目录下")

# 样本信息解析（支持单个/多个样本，从sample_info.txt读取，生成绝对路径）
samples_df = pd.read_csv(
    config["sample_info"],
    sep="\t",
    skiprows=1,
    names=["sample", "raw_fastq"],
    dtype=str
).dropna()
SAMPLES = samples_df["sample"].tolist()
SAMPLE_RAW_FASTQ = dict(zip(samples_df["sample"], samples_df["raw_fastq"]))

# --------------------------
# 最终目标规则（定义所有需生成的最终文件）
# --------------------------
rule all:
    input:
        # 1. GLORI核心输出（m6A位点注释结果）
        GLORI_OUTPUT / f"{SAMPLES[0]}_GLORI_m6A_sites.txt",  # 单个样本简化；多个样本用列表推导
        # 2. 关键中间文件（确保注释和索引构建完成）
        ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base",
        Ref / "hg38.fa",  # 基因组文件
        ANNO_DIR / "GCF_000001405.39_GRCh38.p13_rna2.fa",  # 转录组文件
        # 3. 修剪+去重后的最终Fastq
        TEST_DIR / f"{SAMPLES[0]}_rmdup_trimmed.fq.gz"
    log:
        LOG_DIR / "rule_all_summary.log"


# --------------------------
# 规则1：第一次Trim Galore（原始Fastq修剪）
# --------------------------
rule trim_galore_first:
    input:
        raw_fastq = lambda wildcards: SAMPLE_RAW_FASTQ[wildcards.sample]
    output:
        trimmed_fastq_gz = TEST_DIR / "{sample}_trimmed.fq.gz",
        trim_report = TEST_DIR / "{sample}.fastq.gz_trimming_report.txt"
    params:
        # 从 config.yaml 的 trim_galore_first 节点读取参数
        quality = config["trim_galore_first"]["quality"],
        stringency = config["trim_galore_first"]["stringency"],
        error_rate = config["trim_galore_first"]["error_rate"],
        min_length = config["trim_galore_first"]["min_length"]
    log:
        LOG_DIR / "{sample}_trim_galore_first.log"
    container: config["sif"]
    shell:
        """
            trim_galore -q {params.quality} --stringency {params.stringency} -e {params.error_rate} \
            --length {params.min_length} --fastqc -o {TEST_DIR}  {input.raw_fastq} \
            >> {log} 2>&1
        
        """

# --------------------------
# 规则2：SeqKit去重
# --------------------------
rule seqkit_rmdup:
    input:
        trimmed_fastq_gz = TEST_DIR / "{sample}_trimmed.fq.gz"
    output:
        rmdup_fastq_gz = TEST_DIR / "{sample}_rmdup.fq.gz",
        dup_fq = TEST_DIR / "{sample}_dup.fq.gz"  # 恢复输出，避免变量未定义
    params:
        threads = config["threads"]
    log:
        LOG_DIR / "{sample}_seqkit_rmdup.log"
    container: config["sif"]
    shell:
        """
            seqkit rmdup -j {params.threads} -s -D {output.dup_fq} {input.trimmed_fastq_gz} | gzip > {output.rmdup_fastq_gz} 
        """


# --------------------------
# 规则3：第二次Trim Galore（去重后再次修剪）
# --------------------------
rule trim_galore_second:
    input:
        rmdup_fastq_gz = TEST_DIR / "{sample}_rmdup.fq.gz"
    output:
        final_fastq = TEST_DIR / "{sample}_rmdup_trimmed.fq.gz",
        #final_trim_report = TEST_DIR / "{sample}_rmdup_trimmed_fq.gz_trimming_report.txt"  # 恢复输出，匹配FastQC
    params:
        clip_r1 = config["trim_galore_second"]["clip_r1"],
        quality = config["trim_galore_second"]["quality"],
        min_length = config["trim_galore_second"]["min_length"]
    log:
        LOG_DIR / "{sample}_trim_galore_second.log"
    container: config["sif"]
    shell:
        """
          trim_galore --clip_R1 {params.clip_r1} --quality {params.quality} \
            --length {params.min_length} --fastqc -o {TEST_DIR} {input.rmdup_fastq_gz} \
            >> {log} 2>&1

        """

# --------------------------
# 规则4：转换UCSC GTF格式（change_UCSCgtf.py）
# --------------------------
rule process_ucsc_gtf:
    input:
        gtf = Ref / "GCF_000001405.39_GRCh38.p13_genomic.gtf",  # 依赖提前下载的GTF
        assembly_report = Ref / "GCF_000001405.39_GRCh38.p13_assembly_report.txt"  # 需手动放在ANNO_DIR
    output:
        gtf_converted = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens"  # 直接输出到ANNO_DIR，避免移动
    log:
        LOG_DIR / "process_ucsc_gtf.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/get_anno/change_UCSCgtf.py \
            -i {input.gtf} -j {input.assembly_report} -o {output.gtf_converted} \
            >> {log} 2>&1

        """
# --------------------------
# 规则6：GTF转注释表（gtf2anno.py）
# --------------------------
rule gtf2anno:
    input:
        gtf_converted = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens" 
    output:
        anno_tbl = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl"
    log:
        LOG_DIR / "gtf2anno.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/get_anno/gtf2anno.py \
            -i {input.gtf_converted} -o {output.anno_tbl} \
            >> {log} 2>&1
        
        """

# --------------------------
# 规则7：过滤注释表（awk过滤无效行）
# --------------------------
rule filter_anno_tbl:
    input:
        anno_tbl = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl"  # 依赖规则6输出
    output:
        filtered_tbl = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2"
    log:
        LOG_DIR / "filter_anno_tbl.log"
    container: config["sif"]
    shell:
        """
            awk '$3!~/_/&&$3!="na"' {input.anno_tbl} | sed '/unknown_transcript/d' > {output.filtered_tbl} \
            >> {log} 2>&1
        
        """


# --------------------------
# 规则8：筛选最长转录本（selected_longest_transcrpts_fa.py）
# --------------------------
rule select_longest_transcripts:
    input:
        filtered_tbl = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2",
        rna_fna = Ref / "GCF_000001405.39_GRCh38.p13_rna.fna"  
    output:
        rna_fa = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_rna2.fa"
    log:
        LOG_DIR / "select_longest_transcripts.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/get_anno/selected_longest_transcrpts_fa.py \
            -anno {input.filtered_tbl} -fafile {input.rna_fna} --outname_prx {output.rna_fa} \
            >> {log} 2>&1
        
        """

# --------------------------
# 规则5：构建基因组索引（build_genome_index.py）
# --------------------------
rule build_genome_index:
    input:
        genome_fa = Ref / "hg38.fa"  # 依赖提前下载的基因组
    output:
        genome_index =  INPUT_DIR / "hg38.AG_conversion.fa",  # 索引核心文件
        rvs_genome =  INPUT_DIR / "hg38.rvsCom.fa"  # 反向互补基因组（GLORI需）
    params:
        threads = config["threads"],
        prefix = "hg38"
    log:
        LOG_DIR / "build_genome_index.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/pipelines/build_genome_index.py \
            -f {input.genome_fa} -p {params.threads} -pre {params.prefix} \
            >> {log} 2>&1
 
        """

# --------------------------
# 规则9：构建转录组索引（build_transcriptome_index.py）
# --------------------------
rule build_transcriptome_index:
    input:
        rna_fa = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_rna2.fa"  # 依赖规则8输出
    output:
        transcriptome_index = INPUT_DIR / "GCF_000001405.39_GRCh38.p13_rna2.AG_conversion.fa"
    params:
        threads = config["threads"],
        prefix =  "GCF_000001405.39_GRCh38.p13_rna2"
    log:
        LOG_DIR / "build_transcriptome_index.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/pipelines/build_transcriptome_index.py \
            -p {params.threads} -f {input.rna_fa} -pre {params.prefix} \
            >> {log} 2>&1
        
        """

# --------------------------
# 规则10：生成碱基注释（anno_to_base.py）
# --------------------------
rule anno_to_base:
    input:
        filtered_tbl = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2",
        gtf_converted = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens"  
    output:
        base_anno = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno"
    log:
        LOG_DIR / "anno_to_base.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/get_anno/anno_to_base.py \
            -i {input.filtered_tbl} -o {output.base_anno} \
            >> {log} 2>&1
        
        """

# --------------------------
# 规则11：生成基因列表（gtf2genelist.py）
# --------------------------
rule gtf2genelist:
    input:
        gtf_converted = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens",  # 依赖规则4输出
        rna_fna = Ref / "GCF_000001405.39_GRCh38.p13_rna.fna"  # 提前下载的RNA序列
    output:
        genelist = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist",
        genelist_log = ANNO_DIR / "gtf2genelist_output.log"
    log:
        LOG_DIR / "gtf2genelist.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/get_anno/gtf2genelist.py \
            -i {input.gtf_converted} -f {input.rna_fna} -o {output.genelist} > {output.genelist_log} \
            >> {log} 2>&1
        
        """

# --------------------------
# 规则12：过滤基因列表（awk过滤无效行）
# --------------------------
rule filter_genelist:
    input:
        genelist = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist"  # 依赖规则11输出
    output:
        filtered_genelist = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2"
    log:
        LOG_DIR / "filter_genelist.log"
    container: config["sif"]
    shell:
        """
            awk '$6!~/_/&&$6!="na"' {input.genelist} > {output.filtered_genelist} \
            >> {log} 2>&1
        
        """


# --------------------------
# 规则13：去除注释冗余（anno_to_base_remove_redundance_v1.0.py）
# --------------------------
rule remove_anno_redundance:
    input:
        base_anno = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.baseanno", # 依赖规则10输出
        filtered_genelist = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.genelist2"  # 依赖规则12输出
    output:
        final_base_anno = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base"
    log:
        LOG_DIR / "remove_anno_redundance.log"
    container: config["sif"]
    shell:
        """
            python {TOOL_DIR}/get_anno/anno_to_base_remove_redundance_v1.0.py \
            -i {input.base_anno} -o {output.final_base_anno} -g {input.filtered_genelist} \
            >> {log} 2>&1
        
        """


# --------------------------
# 规则14：运行GLORI调用m6A位点（newrun_GLORI.py）
# --------------------------
rule run_glori:
    input:
        # 输入Fastq（依赖规则3输出）
        final_fastq = TEST_DIR / "{sample}_rmdup_trimmed.fq.gz",
        # 参考文件+索引（依赖规则5、9输出）
        genome_conversion =  INPUT_DIR / "hg38.AG_conversion.fa",
        genome_fa = Ref / "hg38.fa",  # 提前下载的基因组
        rvs_genome =  INPUT_DIR / "hg38.rvsCom.fa",  # 依赖规则5输出
        transcriptome_conversion =  INPUT_DIR / "GCF_000001405.39_GRCh38.p13_rna2.AG_conversion.fa",  # 依赖规则9输出
        # 注释文件（依赖规则7、13输出）
        anno_tbl = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2",
        final_base_anno = ANNO_DIR / "GCF_000001405.39_GRCh38.p13_genomic.gtf_change2Ens.tbl2.noredundance.base"
    output:
        glori_output = GLORI_OUTPUT / "{sample}_GLORI_m6A_sites.txt",
        glori_flag = GLORI_OUTPUT / "{sample}_glori_completed.flag"  # 标记运行完成
        
    params:
        threads = config["threads"],
        tool_dir = TOOL_DIR,
        prefix = "testGLORI_prefix"
    log:
        LOG_DIR / "{sample}_run_glori.log"
    container: config["sif"]
    shell:
        """
            python {params.tool_dir}/newrun_GLORI.py \
            -i {params.tool_dir} \
            -q {input.final_fastq} \
            -T {params.threads} \
            -f {input.genome_conversion} \
            -f2 {input.genome_fa} \
            -rvs {input.rvs_genome} \
            -Tf {input.transcriptome_conversion} \
            -a {input.anno_tbl} \
            -b {input.final_base_anno} \
            -pre {params.prefix} \
            -o {GLORI_OUTPUT} \
            --combine \
            --rvs_fac \
            -c 1 -C 0 -r 0 -p 1.1 -adp 1.1 -s 0 
    
        """