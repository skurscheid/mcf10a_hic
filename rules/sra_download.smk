# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

rule prefetch_se:
    version: 
        1
    conda:
        "../envs/sra_tools.yaml"
    threads:
        1
    log:
        "logs/prefetch/{cell_line}/{chip_antibody}/se/{run}.log"
    input:
    output:
        sra_file = temp("sra_download/{cell_line}/{chip_antibody}/se/{run}.sra")
    shell:
        """
            prefetch {wildcards.run} --output-file {output.sra_file}
        """

rule fastq_dump_se:
    version:
        1
    conda:
        "../envs/sra_tools.yaml"
    threads:
        4
    log:
        "logs/fastq-dump/{cell_line}/{chip_antibody}/se/{run}.log"
    input:
        rules.prefetch_se.output.sra_file
    output:
        temp("raw/{cell_line}/{chip_antibody}/se/{run}.fastq")
    shell:
        """
            fastq-dump --skip-technical {input} --stdout > {output} 2>{log}
        """

rule pigz_fastq_se:
    version:
        1
    threads:
        4
    log:
        "logs/pigz_fastq/{cell_line}/{chip_antibody}/se/{run}.log"
    input:
        rules.fastq_dump_se.output
    output:
        "raw/{cell_line}/{chip_antibody}/se/{run}.fastq.gz"
    shell:
        """
            pigz --processes {threads} --stdout {input} > {output} 2>{log}
        """
    
