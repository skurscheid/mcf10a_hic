__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2018-09-15"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

def get_fastq(wildcards):
    """ returns fastq files for given sample """
    return units.loc[(wildcards["sample"], wildcards["batch"], wildcards["lane"], int(wildcards["replicate"])), ["fq1", "fq2"]].dropna()

def get_fastq_se(wildcards):
    """ returns fastq files for given sample """
    if (wildcards["end"] == "end1"):
      return units.loc[(wildcards["sample"], wildcards["batch"], wildcards["lane"], int(wildcards["replicate"])), ["fq1"]].dropna()
    else:
      return units.loc[(wildcards["sample"], wildcards["batch"], wildcards["lane"], int(wildcards["replicate"])), ["fq2"]].dropna()

singularity: "docker://skurscheid/snakemake_baseimage:0.2"

rule run_fastp_pe:
    conda:
        "../envs/fastqProcessing.yaml"
    version:
        "2"
    threads:
        4
    input:
        get_fastq
    output:
        trimmed_read1 = "fastp/trimmed/pe/{batch}/{sample}_{lane}_{replicate}.end1.fastq.gz",
        trimmed_read2 = "fastp/trimmed/pe/{batch}/{sample}_{lane}_{replicate}.end2.fastq.gz",
        report_html = "fastp/report/pe/{batch}/{sample}_{lane}_{replicate}.fastp.html",
        report_json = "fastp/report/pe/{batch}/{sample}_{lane}_{replicate}.fastp.json"
    shell:
        "fastp -i {input[0]} -I {input[1]} -o {output.trimmed_read1} -O {output.trimmed_read2} --html {output.report_html} --json {output.report_json} --thread {threads}"

rule run_fastp_se:
    conda:
        "../envs/fastqProcessing.yaml"
    version:
        "2"
    threads:
        1
    input:
        get_fastq_se
    output:
        trimmed = "fastp/trimmed/se/{batch}/{sample}_{lane}_{replicate}.{end}.fastq.gz",
        report_html = "fastp/report/se/{batch}/{sample}_{lane}_{replicate}.{end}.fastp.html",
        report_json = "fastp/report/se/{batch}/{sample}_{lane}_{replicate}.{end}.fastp.json"
    shell:
        "fastp -i {input[0]} -o {output.trimmed} --html {output.report_html} --json {output.report_json} --thread {threads}"
