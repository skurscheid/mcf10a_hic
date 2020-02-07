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
        log = "logs/bowtie2_global/{biosample}/{replicate}/{run}_{end}.log"
    input:
        fq = "fastp/trimmed/se/{biosample}/{replicate}/{run}_{end}.fastq.gz"
    output:
        bam = "bowtie2/align_global/se/{biosample}/{replicate}/{run}_{end}.bam",
        unmapped = "bowtie2/align_global/se/{biosample}/{replicate}/{run}_{end}.unmap.fastq"
    shell:
        """
            unset PERL5LIB; bowtie2\
                    -x {params.index}\
                    -p {threads}\
                    -U {input.fq}\
                    {params.cli_params_global}\
                    --un {output.unmapped}\
                    --rg-id BMG\
                    --rg SM:{wildcards.biosample}:{wildcards.run}\
                    2>> {log.log}\
            | samtools view -Shb - > {output.bam}
        """

rule cutsite_trimming:
    """trims potentially chimeric reads prior to second alignment"""
    version:
        1
    params:
        hicpro_dir = config['params']['hicpro']['install_dir']['gadi'],
        cutsite = "AAGCTT" #HindIII
    log:
        log = "logs/cutsite_trimming/{biosample}/{replicate}/{run}_{end}.log"
    input:
        rules.bowtie2_se_global.output.unmapped
    output:
        cutsite_trimmed = "cutsite_trimming/{biosample}/{replicate}/{run}_{end}.fastq"
    shell:
        """ 
            {params.hicpro_dir}/scripts/cutsite_trimming --fastq {input} --cutsite {params.cutsite} --out {output}
        """

rule bowtie2_se_local:
    """ runs alignment of single-end fastq file, modified parameters specific for HiC data"""
    version:
        1
    conda:
        "../envs/fastqProcessing.yaml"
    threads:
        8
    params:
        index = get_index("gadi", config),
        cli_params_local = config['params']['bowtie2']['cli_params_local']
    log:
        log = "logs/bowtie2_local/{biosample}/{replicate}/{run}_{end}.log"
    input:
        fq = rules.cutsite_trimming.output
    output:
        bam = "bowtie2/align_local/se/{biosample}/{replicate}/{run}_{end}.bam",
        unmapped = "bowtie2/align_local/se/{biosample}/{replicate}/{run}_{end}.unmap.fastq"
    shell:
        """
            unset PERL5LIB; bowtie2\
                    -x {params.index}\
                    -p {threads}\
                    -U {input.fq}\
                    {params.cli_params_local}\
                    --un {output.unmapped}\
                    --rg-id BML\
                    --rg SM:{wildcards.biosample}:{wildcards.run}\
                    2>> {log.log}\
            | samtools view -Shb - > {output.bam}
        """
