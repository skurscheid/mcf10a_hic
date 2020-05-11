__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

rule run_fastp_pe:
    conda:
        "../envs/fastqProcessing.yaml"
    version:
        2
    threads:
        4
    log:
        log = "logs/fastp/pe/{biosample}/{replicate}/{run}.log"
    input:
        fq1 = "raw/pe/{cell_line}/{chip_antibody}/{run}_1.fastq.gz",
        fq2 = "raw/pe/{cell_line}/{chip_antibody}/{run}{suffix}.fastq.gz"
    output:
        out1 = "fastp/trimmed/pe/{biosample}/{replicate}/{run}_1.fastq.gz",
        out2 = "fastp/trimmed/pe/{biosample}/{replicate}/{run}_2.fastq.gz",
        report_html = "fastp/report/pe/{biosample}/{replicate}/{run}.fastp.html",
        report_json = "fastp/report/pe/{biosample}/{replicate}/{run}.fastp.json"
    shell:
        """
            fastp -i {input.fq1} -I {input.fq2}\
                  -o {output.out1} -O {output.ou2}\
                  --html {output.report_html} --json {output.report_json}\
                  --length_required 30\
                  --disable_quality_filtering\
                  --detect_adapter_for_pe\
                  --thread {threads} 2>>{log.log}
        """

#rule fastp_dummy:
#    conda:
#        "../envs/fastqProcessing.yaml"
#    version:
#        "2"
#    threads:
#        1
#    input:
#        fastq = "raw/{run}{end}.fastq.gz"
#    output:
#        ln_target = "fastp/trimmed/pe/{biosample}/{replicate}/{run}{end}.fastq.gz"
#    shell:
#        """
#            if [ -e {output.ln_target} ] && [ ! -L {output.ln_target} ];\
#                then rm {output.ln_target}; ln -sr {input.fastq} {output.ln_target};\
#            else\
#                ln -sr {input.fastq} {output.ln_target};\
#            fi
#        """

