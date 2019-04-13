__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2019-04-10"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing HiC data with HiC Explorer
(https://github.com/deeptools/HiCExplorer)

For usage, include this in your workflow.
"""

singularity: "docker://skurscheid/snakemake_baseimage:0.2"

rule findRestSites:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        "1"
    params:
        searchPattern = "AAGCTT"
    input:
        fasta = lambda wildcards: config["params"]["hicexplorer"]["genome_fasta"]["gdu"]
    output:
        rest_sites_bed = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.bed"
    shell:
        """
        hicexplorer --fasta {input.fasta} --searchPattern {params.searchPattern} --outFile {output.rest_sites_bed}
        """

rule hicBuildMatrix_restrictionCutFile_test_run:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
        inputBufferSize = 100000
    threads:
        4
    input:
        mate1 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end1.bam",
        mate2 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end2.bam",
        restrictionCutFile = "hicexplorer/findRestSite/hg38_{resolution}_rest_sites.bed"
    benchmark:
        "hicexplorer/hicBuildMatrix/test_run/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/benchmark/times.tsv"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix/test_run/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix/test_run/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/qc")
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --restrictionCutFile {input.restrictionCutFile} \
                --threads {threads} \
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --QCfolder {output.qcFolder} \
                --doTestRun
        """

rule hicBuildMatrix_restrictionCutFile:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
        inputBufferSize = 100000
    threads:
        4
    input:
        mate1 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end1.bam",
        mate2 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end2.bam",
        restrictionCutFile = "hicexplorer/findRestSite/hg38_{resolution}_rest_sites.bed"
    benchmark:
        "hicexplorer/hicBuildMatrix/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/benchmark/times.tsv"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/qc")
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --restrictionCutFile {input.restrictionCutFile} \
                --threads {threads} \
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --QCfolder {output.qcFolder}
        """
