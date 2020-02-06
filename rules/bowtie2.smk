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

rule bowtie2_se_global:
    """ runs alignment of single-end fastq file, modified parameters specific for HiC data"""
    conda:
        "../envs/fastqProcessing.yaml"
    threads:
        8
    params:
        index = get_index("gadi", config),
        cli_params_global = config['params']['bowtie2']['cli_params_global']
    log:
        log = "logs/bowtie2_global/{sample}_{end}.log"
    input:
        fq = "fastp/trimmed/se/{sample}_{end}.fastq.gz"
    output:
        bam = "bowtie2/align_global/se/{sample}_{end}.bam",
        unmapped = "bowtie2/align_global/se/{sample}_{end}.unmap.fastq"
    shell:
        """
            unset PERL5LIB; bowtie2\
                    -x {params.index}\
                    -p {threads}\
                    -U {input.fq}\
                    {params.cli_params}\
                    --un {output.unmapped}\
                    -rg-id BMG\
                    --rg SM:{wildcards.sample}\
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
        index = get_index("gadi", config),
        cli_params = config['params']['bowtie2']['cli_params_local']
    log:
        log = "logs/bowtie2_local//{sample}_{end}.log"
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
