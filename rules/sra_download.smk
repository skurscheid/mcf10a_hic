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
            prefetch {wildcards.run} --output-file {output.sra_file} 2>{log}
        """

rule prefetch_pe:
    version:
        1
    conda:
        "../envs/sra_tools.yaml"
    threads:
        1
    params:
        max_size = lambda wildcards: list(runTable[runTable.Run == wildcards['run']].Bytes + 1000000)
    log:
        "logs/prefetch/pe/{cell_line}/{chip_antibody}/{run}.log"
    input:
    output:
        sra_file = temp("sra_download/pe/{cell_line}/{chip_antibody}/{run}.sra")
    shell:
        """
            prefetch --max-size {params.max_size} {wildcards.run} --output-file {output.sra_file} 2>{log}
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

rule fastq_dump_pe:
    version:
        1
    conda:
        "../envs/sra_tools.yaml"
    threads:
        4
    log:
        "logs/fastq-dump/{cell_line}/{chip_antibody}/pe/{run}.log"
    input:
        rules.prefetch_pe.output.sra_file
    output:
        temp(directory("raw_fastq/pe/{cell_line}/{chip_antibody}/{run}")),
        "raw_fastq/pe/{cell_line}/{chip_antibody}/{run}/{run}_1.fastq",
        "raw_fastq/pe/{cell_line}/{chip_antibody}/{run}/{run}_2.fastq"
    shell:
        """
            fastq-dump --skip-technical --split-files {input} --outdir {output[0]} 2>{log}
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
    
rule pigz_fastq_pe:
    version:
        1
    threads:
        4
    log:
        "logs/pigz_fastq/pe/{cell_line}/{chip_antibody}/{run}{suffix}.log"
    input:
        dir = rules.fastq_dump_pe.output,
        file = lambda wildcards: "/".join(['raw_fastq', 'pe', wildcards['cell_line'], wildcards['chip_antibody'], wildcards['run'], wildcards['run'] + wildcards['suffix'] + '.fastq'])
    output:
        "raw/pe/{cell_line}/{chip_antibody}/{run}/{run}{suffix}.fastq.gz"
    shell:
        """
            pigz --processes {threads} --stdout {input.file} > {output} 2>{log}
        """
    
