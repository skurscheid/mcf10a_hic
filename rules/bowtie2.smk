__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2019-04-08"

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
        index = get_index("gdu", config),
        cli_params = config['params']['bowtie2']['cli_params']
    input:
        fq = "fastp/trimmed/se/{batch}/{sample}_{lane}_{replicate}.{end}.fastq.gz"
    output:
        bam = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.{end}.bam",
        metrics = "bowtie2/report/se/{batch}/{sample}_{lane}_{replicate}.{end}.txt"
    shell:
        """
            unset PERL5LIB; bowtie2 {params.cli_params}\
                    -p {threads}\
                    -x {params.index}\
                    -U {input.fq}\
                    --rg-id {wildcards.sample}:{wildcards.batch}\
                    --rg "replicate:{wildcards.replicate}"\
                    --met-file {output.metrics}\
            | samtools view -Sb - > {output.bam}
        """
