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
        searchPattern = lambda wildcards: config["params"]["hicexplorer"]["findRestSites"][wildcards.res_enzyme]
    input:
        fasta = lambda wildcards: config["params"]["hicexplorer"]["genome_fasta"]["gadi"]
    output:
        rest_sites_bed = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.bed"
    shell:
        """
            findRestSite --fasta {input.fasta} --searchPattern {params.searchPattern} --outFile {output.rest_sites_bed}
        """

rule mappableRestSites:
    version:
        "1"
    params:
    input:
        rest_sites_bed = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.bed"
    output:
        mappable_rest_sites_bed = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.{kmer}.bed"
    shell:
        "touch {output.mappable_rest_sites_bed}"

rule hicBuildMatrix_restrictionCutFile_test_run:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        5
    params:
        inputBufferSize = 400000
    threads:
        8
    input:
        mate1 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end1_suffix"] + ".bam",
        mate2 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end2_suffix"] + ".bam",
        restrictionCutFile = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.bed"
    benchmark:
        "benchmarks/hicexplorer/hicBuildMatrix_rest/test_run/{res_enzyme}/{biosample}/{rep}/{run}/times.tsv"
    log: 
        "logs/hicexplorer/hicBuildMatrix_rest/test_run/{res_enzyme}/{biosample}/{rep}/{run}/log.txt"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_rest/test_run/{res_enzyme}/{biosample}/{rep}/{run}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix_rest/test_run/{res_enzyme}/{biosample}/{rep}/{run}/qc"),
        outBam = "hicexplorer/hicBuildMatrix_rest/test_run/{res_enzyme}/{biosample}/{rep}/{run}_hic_matrix.bam"
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --restrictionCutFile {input.restrictionCutFile} \
                --threads {threads} \
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --outBam {output.outBam}\
                --QCfolder {output.qcFolder} \
                --doTestRun 1>{log} 2>{log}
        """

rule hicBuildMatrix_restrictionCutFile:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        5
    params:
        inputBufferSize = 400000
    threads:
        8
    input:
        mate1 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end1_suffix"] + ".bam",
        mate2 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end2_suffix"] + ".bam",
        restrictionCutFile = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.bed"
    benchmark:
        "benchmarks/hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{biosample}/{rep}/{run}/times.tsv"
    log: 
        "logs/hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{biosample}/{rep}/{run}/log.txt"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{biosample}/{rep}/{run}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{biosample}/{rep}/{run}/qc"),
        outBam = "hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{biosample}/{rep}/{run}_hic_matrix.bam"
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --restrictionCutFile {input.restrictionCutFile} \
                --threads {threads} \
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --outBam {output.outBam}\
                --QCfolder {output.qcFolder} \
                1>{log} 2>{log}
        """

rule hicBuildMatrix_bin_test_run:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
        inputBufferSize = 400000,
        restrictionSequence = "AAGCTT"
    threads:
        8
    input:
        mate1 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end1_suffix"] + ".bam",
        mate2 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end2_suffix"] + ".bam"
    benchmark:
        "benchmarks/hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}/times.tsv"
    log:
        "logs/hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}/log.txt"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}/qc"),
        outBam = "hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}_hic_matrix.bam"
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --threads {threads}\
                --restrictionSequence {params.restrictionSequence}\
                --binSize {wildcards.resolution}\
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --outBam {output.outBam}\
                --QCfolder {output.qcFolder} \
                --doTestRun 1>{log} 2>{log}
        """

rule hicBuildMatrix_bin:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        4
    params:
        inputBufferSize = 400000
    threads:
        8
    input:
        mate1 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end1_suffix"] + ".bam",
        mate2 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end2_suffix"] + ".bam"
    log:
        "logs/hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}.txt"
    benchmark:
        "benchmarks/hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}.tsv"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix_bin/{resolution}/{biosample}/{rep}/{run}/qc")
    shell:
        """
            hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                    --threads {threads} \
                    --binSize {wildcards.resolution}\
                    --inputBufferSize {params.inputBufferSize} \
                    --outFileName {output.outHicMatrix} \
                    --QCfolder {output.qcFolder} 1>{log} 2>{log}
        """

