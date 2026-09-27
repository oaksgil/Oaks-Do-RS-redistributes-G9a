#########################################
# Snakemake pipeline for C&R analysis
# Brian Do
# 220612
# Adapted from https://github.com/BleekerLab/snakemake_rnaseq 
#########################################


###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############
configfile: "config/bd523_bd524_cutandrun.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]


########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information
samples = pd.read_csv(config["sample_sheet"], dtype=str).set_index(["name"], drop=False)
SAMPLES = samples.index.tolist()

###########################
# Input functions for rules
###########################

def get_fastq(wildcards):
    return({"forward_read1": config["fastq_dir"] + samples.loc[(wildcards.sample), "filename_r1s1"],
            "reverse_read1": config["fastq_dir"] + samples.loc[(wildcards.sample), "filename_r2s1"],
            "forward_read2": config["fastq_dir"] + samples.loc[(wildcards.sample), "filename_r1s2"],
            "reverse_read2": config["fastq_dir"] + samples.loc[(wildcards.sample), "filename_r2s2"]})

def extract_mapped_reads(bowtie_file_path):
    """
    Extract the sum of mapped reads from a bowtie2 log file.
    
    Args:
        bowtie_file_path (str): Path to the bowtie2 log file
        
    Returns:
        int: Sum of mapped reads (concordantly aligned exactly 1 time + >1 times)
        float: Overall alignment rate percentage
    """
    with open(bowtie_file_path, 'r') as file:
        content = file.read()
        
        # Extract reads aligned exactly 1 time
        aligned_once_line = [line for line in content.split('\n') if "aligned concordantly exactly 1 time" in line][0]
        aligned_once = int(aligned_once_line.split()[0])
        
        # Extract reads aligned more than 1 time
        aligned_multiple_line = [line for line in content.split('\n') if "aligned concordantly >1 times" in line][0]
        aligned_multiple = int(aligned_multiple_line.split()[0])
        
        # Calculate sum of mapped reads
        sum_mapped_reads = aligned_once + aligned_multiple

        return sum_mapped_reads


#################
# Desired outputs
#################
TRIMMED_FILES = expand(RESULT_DIR + "fastq_trimmed/{sample}.read1.fastq.gz", sample = SAMPLES)
BAM_FILES = expand(RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup_noblacklist.bam", sample = SAMPLES)
MOUSE_BAM = expand(RESULT_DIR + "bowtie2_bam/{sample}_mouse_aligned_reads.bam", sample = SAMPLES)
BW_FILES = expand(RESULT_DIR + "bigwigs/{sample}_norm.bw", sample = SAMPLES)
NO_RMDUP_BW_FILES = expand(RESULT_DIR + "bigwigs/{sample}_no_rmdup_norm.bw", sample = SAMPLES)
NO_RMDUP_MS_BW_FILES = expand(RESULT_DIR + "bigwigs/norm_to_mouse/{sample}_no_rmdup_norm_to_mouse.bw", sample = SAMPLES)
FASTQC_FILES1 = expand(RESULT_DIR + "fastq_trimmed/fastqc/{sample}.read1_fastqc.zip", sample=SAMPLES)
FASTQC_FILES2 = expand(RESULT_DIR + "fastq_trimmed/fastqc/{sample}.read2_fastqc.zip", sample=SAMPLES)

rule all:
    input:
        TRIMMED_FILES,
        BAM_FILES, 
        MOUSE_BAM,
        BW_FILES,
        NO_RMDUP_BW_FILES,
        NO_RMDUP_MS_BW_FILES,
        "multiqc_report.html"
    message:
        "CUT&RUN run complete!"

#######
# Rules
#######

# Removing adapter sequences
# Sometimes read 2 is entirely G's
# Sometimes the entire read is adapter (for some reason without an A)

rule trim_filter:
    input:
        unpack(get_fastq)
    output:
        fwd = RESULT_DIR + "fastq_trimmed/{sample}.read1.fastq.gz",
        rev = RESULT_DIR + "fastq_trimmed/{sample}.read2.fastq.gz"
    params:
        sample_name =  "{sample}",
        threads     =  config["threads"],
        working_dir = WORKING_DIR,
        tmp = WORKING_DIR + "{sample}"
    shell: 
        "mkdir -p {params.working_dir}; \
        cat {input.forward_read1} {input.forward_read2} > {params.tmp}_fwd.fastq.gz; \
        cat {input.reverse_read1} {input.reverse_read2} > {params.tmp}_rev.fastq.gz; \
        cutadapt -j {params.threads} --nextseq-trim=20 --minimum-length 20 \
        -a AGATCGGAAGAGCACA -A AGATCGGAAGAGCGTC \
        -g ^GATCGGAAGAGCACA -G ^GATCGGAAGAGCG -G GGGGGGGGGGG \
        -o {output.fwd} -p {output.rev} \
         {params.tmp}_fwd.fastq.gz {params.tmp}_rev.fastq.gz; \
        rm -f {params.tmp}_fwd.fastq.gz {params.tmp}_rev.fastq.gz"

# Map to human then mouse

rule map_to_human_genome_using_bowtie:
    input:
        fwd = RESULT_DIR + "fastq_trimmed/{sample}.read1.fastq.gz",
        rev = RESULT_DIR + "fastq_trimmed/{sample}.read2.fastq.gz"
    output:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads.bam",
        RESULT_DIR + "bowtie2_bam/{sample}_unaligned.fastq.1.gz",
        RESULT_DIR + "bowtie2_bam/{sample}_unaligned.fastq.2.gz"
    message:
        "mapping {wildcards.sample} reads to human genome"
    params:
        sample_name           =  "{sample}",
        genome_index          =  config["hs_genome_index"],
        prefix                =  RESULT_DIR + "bowtie2_bam/{sample}",
        threads               =  config["threads"]
    threads: int(config["threads"])
    resources: 
        cpus = config["threads"]
    shell:
        "mkdir -p logs; "
        "bowtie2 -p {params.threads} --dovetail -I 10 -X 700 \
        --no-mixed --no-discordant --very-sensitive-local --local \
        --phred33 -x {params.genome_index} --un-conc-gz {params.prefix}_unaligned.fastq.gz \
        -1 {input.fwd} -2 {input.rev} 2> logs/{wildcards.sample}_log.bowtie2 | \
        samtools view -bS - > {params.prefix}_aligned_reads.bam"

rule map_to_mouse_genome_using_bowtie:
    input:
        fwd = RESULT_DIR + "bowtie2_bam/{sample}_unaligned.fastq.1.gz",
        rev = RESULT_DIR + "bowtie2_bam/{sample}_unaligned.fastq.2.gz"
    output:
        RESULT_DIR + "bowtie2_bam/{sample}_mouse_aligned_reads.bam"
    message:
        "mapping {wildcards.sample} unmapped reads to mouse genome"
    params:
        sample_name           =  "{sample}",
        genome_index          =  config["mm_genome_index"],
        prefix                =  RESULT_DIR + "bowtie2_bam/{sample}",
        threads               =  config["threads"]
    threads: int(config["threads"])
    resources: 
        cpus = config["threads"]
    shell:
        "mkdir -p logs; "
        "bowtie2 -p {params.threads} --dovetail -I 10 -X 700 \
        --no-mixed --no-discordant --very-sensitive-local --local \
        --phred33 --no-unal -x {params.genome_index} -1 {input.fwd} -2 {input.rev} \
        2> logs/{wildcards.sample}_mouse_log.bowtie2 | \
        samtools view -bS - > {params.prefix}_mouse_aligned_reads.bam"

# Duplicate removal

rule sort_rmdup_bowtie_and_index:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads.bam"
    output:
        sorted_output = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup.bam",
        sorted_index = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup.bam.bai"
    params:
        working_dir = WORKING_DIR,
        tmp = WORKING_DIR + "{sample}1"
    shell:
        "mkdir -p {params.working_dir}; samtools view -h -f 0x2 {input} | samtools collate -o - - | samtools fixmate -m - - | \
        samtools sort -T {params.tmp} -o - - | samtools markdup -r - {output.sorted_output} && samtools index {output.sorted_output}"

rule rmdup_remove_blacklist:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup.bam"
    output:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup_noblacklist.bam"
    params:
        blacklist = config["hs_blacklist"]
    shell:
        "bedtools intersect -v -abam {input} -b {params.blacklist} > {output} && samtools index {output}"

rule rmdup_normalize_cpm_and_make_tracks:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup_noblacklist.bam"
    output:
        RESULT_DIR + "bigwigs/{sample}_norm.bw"
    params:
        threads = config["threads"],
        aligned_rmdup = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_rmdup.bam"
    shell:
        "bamCoverage -b {input} -o {output} --normalizeUsing CPM -p {params.threads} -bs 1 && "
        "rm -f {params.aligned_rmdup} && rm -f {params.aligned_rmdup}.bai"

# No duplicate removal

rule no_rmdup_index:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads.bam"
    output:
        sorted_output = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted.bam",
        sorted_index = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted.bam.bai"
    params:
        working_dir = WORKING_DIR,
        tmp = WORKING_DIR + "{sample}2"
    shell:
        "mkdir -p {params.working_dir}; "
        "samtools view -h -f 0x2 {input} | samtools sort -T {params.tmp} -o {output.sorted_output} - && samtools index {output.sorted_output}"

rule no_rmdup_remove_blacklist:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted.bam"
    output:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_noblacklist.bam"
    params:
        blacklist = config["hs_blacklist"]
    shell:
        "bedtools intersect -v -abam {input} -b {params.blacklist} > {output} && samtools index {output}"

rule no_rmdup_normalize_cpm_and_make_tracks:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_noblacklist.bam"
    output:
        RESULT_DIR + "bigwigs/{sample}_no_rmdup_norm.bw"
    params:
        threads = config["threads"],
        aligned_rmdup = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted.bam"
    shell:
        "bamCoverage -b {input} -o {output} --normalizeUsing CPM -p {params.threads} -bs 1 && "
        "rm -f {params.aligned_rmdup} && rm -f {params.aligned_rmdup}.bai"

# Normalize to # of mouse reads

rule no_rmdup_normalize_to_mouse_and_make_tracks:
    input:
        bam = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted_noblacklist.bam",
        human_log = "logs/{sample}_log.bowtie2",
        mouse_log = "logs/{sample}_mouse_log.bowtie2"
    output:
        RESULT_DIR + "bigwigs/norm_to_mouse/{sample}_no_rmdup_norm_to_mouse.bw"
    params:
        threads = config["threads"],
        aligned_rmdup = RESULT_DIR + "bowtie2_bam/{sample}_aligned_sorted.bam",
        output_dir = RESULT_DIR + "bigwigs/norm_to_mouse/",
        eff_genome_sz = 2913022398
    run:
        human_reads = extract_mapped_reads(input.human_log)
        mouse_reads = extract_mapped_reads(input.mouse_log)
        sf = human_reads/mouse_reads/100

        shell_cmd = "mkdir -p {params.output_dir}; "
        shell_cmd += f"bamCoverage -b {input.bam} -o {output} --normalizeUsing RPGC \
        --ignoreForNormalization chrX chrM --effectiveGenomeSize {params.eff_genome_sz} \
        --scaleFactor {sf} -p {params.threads} -bs 10 && "
        shell_cmd += f"rm -f {params.aligned_rmdup} && rm -f {params.aligned_rmdup}.bai"

        print(shell_cmd)
        shell(shell_cmd)

# QC

rule run_fastqc:
    input:
        fwd = RESULT_DIR + "fastq_trimmed/{sample}.read1.fastq.gz",
        rev = RESULT_DIR + "fastq_trimmed/{sample}.read2.fastq.gz",
    output:
        RESULT_DIR + "fastq_trimmed/fastqc/{sample}.read1_fastqc.zip",
        RESULT_DIR + "fastq_trimmed/fastqc/{sample}.read2_fastqc.zip",
    params:
        out_dir = RESULT_DIR + "fastq_trimmed/fastqc/"
    shell:
        "mkdir -p {params.out_dir}; \
         fastqc -o {params.out_dir} {input.fwd} & \
         fastqc -o {params.out_dir} {input.rev}"

rule run_multiqc:
    input:
        BAM_FILES,
        FASTQC_FILES1,
        FASTQC_FILES2
    output:
        multiqc = "multiqc_report.html"
    shell:
        "multiqc . --ignore \"logs/*mouse*\" --ignore \"results/bowtie2_bam/*mouse*\""







