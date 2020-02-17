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

#rule run_fastp_se:
#    conda:
#        "../envs/fastqProcessing.yaml"
#    version:
#        "2"
#    threads:
#        4
#    input:
#        fastq = "raw/{run}{end}.fastq.gz"
#   output:
#        trimmed = "fastp/trimmed/se/{biosample}/{replicate}/{run}{end}.fastq.gz",
#        report_html = "fastp/report/se/{biosample}/{replicate}/{run}{end}.fastp.html",
#        report_json = "fastp/report/se/{biosample}/{replicate}/{run}{end}.fastp.json"
#    shell:
#        "fastp -i {input[0]} -o {output.trimmed} --html {output.report_html} --json {output.report_json} --thread {threads}"

rule fastp_dummy:
    conda:
        "../envs/fastqProcessing.yaml"
    version:
        "2"
    threads:
        1
    input:
        fastq = "raw/{run}{end}.fastq.gz"
    output:
        ln_target = "fastp/trimmed/se/{biosample}/{replicate}/{run}{end}.fastq.gz"
    shell:
        """
            if [ -e {output.ln_target} ] && [ ! -L {output.ln_target} ];\
                then rm {output.ln_target}; ln -sr {input.fastq} {output.ln_target};\
            else\
                ln -sr {input.fastq} {output.ln_target};\
            fi
        """

