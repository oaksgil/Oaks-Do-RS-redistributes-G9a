#########################################
# Snakemake pipeline for ATACseq analysis
# Brian Do
# 210407
# Adapted from https://github.com/BleekerLab/snakemake_rnaseq 
#########################################


###########
# Libraries
###########
import pandas as pd

###############
# Configuration
###############
configfile: "config/220526_bd186_hudep_atacseq.yaml" # where to find parameters
WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]


########################
# Samples and conditions
########################

# read the tabulated separated table containing the sample, condition and fastq file information∂DE
samples = pd.read_csv(config["sample_sheet"], dtype=str).set_index(["name"], drop=False)
SAMPLES = samples.index.tolist()

###########################
# Input functions for rules
###########################

def get_fastq(wildcards):
    return({"forward_read": config["fastq_dir"] + samples.loc[(wildcards.sample), "filename1"],
            "reverse_read": config["fastq_dir"] + samples.loc[(wildcards.sample), "filename2"]})

def get_norm_factor(wildcards, input):
    if not Path(input.norm_factors_file).exists():
        return(-1)

    nf = pd.read_table(input.norm_factors_file)  # see the normalize_bedgraphs rule
    return(nf.loc[wildcards.sample, "div.factors"])

#################
# Desired outputs
#################
FASTQC_FILES = expand("logs/{sample}_r1_trimmed_fastqc.html", sample=SAMPLES)
BAM_FILES = expand(RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads_sorted_nodup_noblacklist_noM.bam", sample = SAMPLES)
HISTO_FILES = expand("logs/{sample}_picard.pdf", sample=SAMPLES)
BDG_FILES = expand(RESULT_DIR + "macs2/{sample}_treat_pileup.bdg", sample = SAMPLES)
NORM_BDG_FILES = expand(RESULT_DIR + "norm_bedgraphs/{sample}_norm.bedgraph", sample=SAMPLES)
BW_FILES = expand(RESULT_DIR + "norm_bigwigs/{sample}.bw", sample = SAMPLES)

rule all:
    input:
        FASTQC_FILES,
        BAM_FILES, 
        HISTO_FILES,
        BDG_FILES,
        RESULT_DIR + "counts/" + config["peak_name"] + "_intersection_peaks.saf",
        RESULT_DIR + "counts/" + config["peak_name"] + "_intersection_counts.txt",
        RESULT_DIR + "counts/" + config["peak_name"] + "_norm_factors.txt",
        NORM_BDG_FILES,
        BW_FILES
    message:
        "ATACseq run complete!"

#######
# Rules
#######

rule trim_adapters:
    input:
        unpack(get_fastq)
    output:
        r1 = WORKING_DIR  + "{sample}_r1_trimmed.fastq.gz",
        r2 = WORKING_DIR  + "{sample}_r2_trimmed.fastq.gz"
    params:
        threads =  config["threads"]
    shell:
        "mkdir -p logs; mkdir -p plots; cutadapt -j {params.threads} -m 20 -q 20 -a CTGTCTCTTATA -A CTGTCTCTTATA -o {output.r1} -p {output.r2} \
        {input.forward_read} {input.reverse_read} > logs/{wildcards.sample}_cutadapt.log"

rule fastqc:
    input:
        r1 = WORKING_DIR  + "{sample}_r1_trimmed.fastq.gz",
        r2 = WORKING_DIR  + "{sample}_r2_trimmed.fastq.gz"
    output:
        "logs/{sample}_r1_trimmed_fastqc.html",
        "logs/{sample}_r2_trimmed_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} -o logs"

rule map_to_genome_using_bowtie:
    input:
        r1 = WORKING_DIR  + "{sample}_r1_trimmed.fastq.gz",
        r2 = WORKING_DIR  + "{sample}_r2_trimmed.fastq.gz"
    output:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads_sorted_nodup_noblacklist_noM.bam"
    message:
        "mapping {wildcards.sample} reads to genome"
    params:
        sample_name           =  "{sample}",
        genome_index          =  config["genome_index"],
        prefix                =  WORKING_DIR + "{sample}",
        threads               =  config["threads"],
        blacklist             =  config["blacklist"],
        tmp                   =  WORKING_DIR + "{sample}"
    threads: int(config["threads"])
    resources: 
        cpus = config["threads"]
    shell:  # aligning end-to-end with max fragment length 1000 and min 10, then fix mates, remove chrM and blacklist, sort, and remove duplicates 
        "mkdir -p logs; "
        "bowtie2 -p {params.threads} --dovetail -I 10 -X 1000 --very-sensitive \
        --phred33 --no-unal -x {params.genome_index} -1 {input.r1} -2 {input.r2} 2> logs/{wildcards.sample}_log.bowtie2 | \
        samtools fixmate -O sam -m - - | \
        grep -v 'chrM\|alt\|random\|chrUn' | \
        samtools view -f 2 -ub - | \
        bedtools intersect -v -ubam -abam - -b {params.blacklist} | \
        samtools sort -u -@{params.threads} -T {params.tmp} - | \
        samtools markdup -r -s -@{params.threads} - {output} 2> logs/{wildcards.sample}_markdup.log"

# picard collect insert size metrics 
rule picard:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads_sorted_nodup_noblacklist_noM.bam"
    output:
        log = "logs/{sample}_picard.txt",
        histo = "logs/{sample}_picard.pdf"
    shell:
        "picard CollectInsertSizeMetrics -I {input} -O {output.log} -H {output.histo} -W 1000"

# MACS2 shifting and peak calling
rule macs2:
    input:
        RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads_sorted_nodup_noblacklist_noM.bam"
    output:
        RESULT_DIR + "macs2/{sample}_peaks.narrowPeak",
        RESULT_DIR + "macs2/{sample}_treat_pileup.bdg",
        RESULT_DIR + "macs2/{sample}_summits.bed"
    params:
        bed = WORKING_DIR + "{sample}.bed",
        genome_size = config["macs2_genome_size"],
        resdir = RESULT_DIR + "macs2/"
    shell:
        "mkdir -p {params.resdir}; bedtools bamtobed -i {input} > {params.bed}; \
        macs2 callpeak -f BED -t {params.bed} -g {params.genome_size} \
        -B --nomodel --shift -100 --extsize 200 --keep-dup all --call-summits -n {wildcards.sample} \
        --outdir {params.resdir} 2> logs/{wildcards.sample}.macs2.log"

# get a list of peaks
rule get_expanded_peaks:  # need to shift and extend
    input:
        RESULT_DIR + "macs2/{sample}_summits.bed"
    output:
        WORKING_DIR + "{sample}_summits_pm500.bed"
    shell:
        "awk 'BEGIN {{OFS=\"\\t\"}}; {{print($1,$2-499,$3+499,$4,$5)}}' {input} | sort -k5,5nr > {output}"

# get a list of peaks that are shared by all files 
rule intersect_peaks: 
    input:
        bdg_files = [WORKING_DIR + "{sample}_summits_pm500.bed".format(sample=sample) for sample in SAMPLES]
    output:
        RESULT_DIR + "counts/" + config["peak_name"] + "_intersection_peaks.saf"
    shell:
        "bedtools intersect -a {input.bdg_files[0]} -b {input.bdg_files} -wa -u | sort -k1,1 -k2,2n | \
        awk 'BEGIN {{OFS=\"\\t\"; print(\"GeneID\", \"Chr\", \"Start\", \"End\", \"Strand\")}}; \
        {{print($4,$1,$2,$3,\"+\")}}' > {output}" 

# count features in the shared SAF file
rule count_features:
    input:
        saf = RESULT_DIR + "counts/" + config["peak_name"] + "_intersection_peaks.saf",
        bams = expand(RESULT_DIR + "bowtie2_bam/{sample}_aligned_reads_sorted_nodup_noblacklist_noM.bam", sample = SAMPLES)
    output:
        RESULT_DIR + "counts/" + config["peak_name"] + "_intersection_counts.txt"
    params:
        threads = config["threads"]
    shell:
        "featureCounts -T {params.threads} -M -O -F 'SAF' -a {input.saf} -o {output} {input.bams}"

# get normalization factors in R with edgeR and TMM
rule get_norm_factors:
    input:
        cts = RESULT_DIR + "counts/" + config["peak_name"] + "_intersection_counts.txt",
        script = "scripts/atacseq_normalization.r"
    output:
        pdf = "plots/" + config["peak_name"] + "_mds.pdf",
        txt = RESULT_DIR + "counts/" + config["peak_name"] + "_norm_factors.txt"
    shell:
        "mkdir -p plots; Rscript {input.script} {input.cts} {output.pdf} {output.txt}"

# normalize bedgraphs
rule normalize_bedgraphs:
    input:
        norm_factors_file = RESULT_DIR + "counts/" + config["peak_name"] + "_norm_factors.txt",
        bdg = RESULT_DIR + "macs2/{sample}_treat_pileup.bdg"
    params:
        div_factor = get_norm_factor  # takes in wildcards and input.norm_factors_file
    output:
        RESULT_DIR + "norm_bedgraphs/{sample}_norm.bedgraph"
    shell:
        "awk 'BEGIN {{OFS=\"\\t\"}}; {{printf(\"%s %s %s %.3f\\n\", \
        $1,$2,$3,$4/{params.div_factor})}}' {input.bdg} > {output}"

# make bigwig 
rule make_bigwigs:
    input:
        RESULT_DIR + "norm_bedgraphs/{sample}_norm.bedgraph"
    output:
        RESULT_DIR + "norm_bigwigs/{sample}.bw"
    params:
        chromsize = config["chromsize"]
    shell:
        "bedGraphToBigWig {input} {params.chromsize} {output}"























