__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for aligning reads with bowtie2
(https://github.com/BenLangmead/bowtie2)

For usage, include this in your workflow.
"""

def get_index(machine, config):
    """ returns path to index"""
    return config["params"]["bowtie2"]["index"][machine]

singularity: "docker://skurscheid/snakemake_baseimage:0.2"

rule bowtie2_se:
    """ runs alignment of single-end fastq file, modified parameters specific for HiC data"""
    conda:
        "../envs/fastqProcessing.yaml"
    threads:
        8
    params:
        index = get_index("gadi", config),
        cli_params = config['params']['bowtie2']['cli_params']
    log:
        log = "logs/bowtie2/{sample}_{end}.log",
        metrics = "bowtie2/report/se/{sample}_{end}.txt"
    input:
        fq = "fastp/trimmed/se/{sample}_{end}.fastq.gz"
    output:
        bam = "bowtie2/align/se/{sample}_{end}.bam"
    shell:
        """
            unset PERL5LIB; bowtie2\
                    -x {params.index}\
                    -p {threads}\
                    -U {input.fq}\
                    {params.cli_params}\
                    --rg-id {wildcards.sample}\
                    --met-file {log.metrics}\
                    2 >> {log.log}\
            | samtools view -Shb - > {output.bam}
        """

rule bowtie2_se_rerun:
    """ runs alignment of single-end fastq file, modified parameters specific for HiC data"""
    version:
        1
    conda:
        "../envs/fastqProcessing.yaml"
    threads:
        8
    params:
        index = get_index("gdu", config),
        cli_params = "--reorder"
    log:
        log = "logs/bowtie2_se_rerun/{sample}_{end}.log",
        metrics = "bowtie2_se_rerun/report/se/{sample}_{end}.txt"
    input:
        fq = "fastp/trimmed/pe/{sample}_{end}.fastq.gz"
    output:
        bam = "bowtie2_rerun/align/se/{sample}_{end}.bam",
        metrics = "bowtie2_rerun/report/se/{sample}_{end}.txt"
    shell:
        """
            unset PERL5LIB; bowtie2\
                    -x {params.index}\
                    -p {threads}\
                    -U {input.fq}\
                    {params.cli_params}\
                    --met-file {log.metrics}\
                    2 >> {log.log}\
            | samtools view -Shb - > {output.bam}
        """
