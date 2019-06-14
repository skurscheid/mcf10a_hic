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
        mate1 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end1.bam",
        mate2 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end2.bam",
        restrictionCutFile = "hicexplorer/findRestSite/hg38_{res_enzyme}_rest_sites.k50.bed"
    benchmark:
        "hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/test_{sample}_{lane}_{replicate}/benchmark/times.tsv"
    log: 
        "hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/test_{sample}_{lane}_{replicate}/log.txt"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/test_{sample}_{lane}_{replicate}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/test_{sample}_{lane}_{replicate}/qc"),
        outBam = "hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/test_{sample}_{lane}_{replicate}_hic_matrix.bam"
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
        mate1 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end1.bam",
        mate2 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end2.bam"
    benchmark:
        "hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{batch}/{sample}/test_{sample}_{lane}_{replicate}/benchmark/times.tsv"
    log:
        "hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{batch}/{sample}/test_{sample}_{lane}_{replicate}/log.txt"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{batch}/{sample}/test_{sample}_{lane}_{replicate}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{batch}/{sample}/test_{sample}_{lane}_{replicate}/qc"),
        outBam = "hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{batch}/{sample}/test_{sample}_{lane}_{replicate}_hic_matrix.bam"
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

rule hicBuildMatrix_restrictionCutFile:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        4
    params:
        inputBufferSize = 400000
    threads:
        8
    input:
        mate1 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end1.bam",
        mate2 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end2.bam",
        restrictionCutFile = "hicexplorer/findRestSite/hg38_{rest_site}_rest_sites.k50.bed"
    benchmark:
        "hicexplorer/hicBuildMatrix/{rest_site}/{batch}/{sample}/{sample}_{lane}_{replicate}/benchmark/times.tsv"
    log:
        "hicexplorer/hicBuildMatrix/{rest_site}/{batch}/{sample}/{sample}_{lane}_{replicate}/log.txt"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix/{rest_site}/{batch}/{sample}/{sample}_{lane}_{replicate}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix/{rest_site}/{batch}/{sample}/{sample}_{lane}_{replicate}/qc")
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --restrictionCutFile {input.restrictionCutFile} \
                --threads {threads} \
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --QCfolder {output.qcFolder} 1>{log} 2>{log}
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
        mate1 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end1.bam",
        mate2 = "bowtie2/align/se/{batch}/{sample}_{lane}_{replicate}.end2.bam"
    log:
        "hicexplorer/hicBuildMatrix_bin/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/log.txt"
    benchmark:
        "hicexplorer/hicBuildMatrix_bin/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/benchmark/times.tsv"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_bin/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}_hic_matrix.h5",
        qcFolder = directory("hicexplorer/hicBuildMatrix_bin/{resolution}/{batch}/{sample}/{sample}_{lane}_{replicate}/qc")
    shell:
        """
        hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                --threads {threads} \
                --binSize {wildcards.resolution}\
                --inputBufferSize {params.inputBufferSize} \
                --outFileName {output.outHicMatrix} \
                --QCfolder {output.qcFolder} 1>{log} 2>{log}
        """

rule hicQC_per_batch:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
        labels = hicQCLabels
    threads:
        64
    input:
        hicQCInput
    output:
        directory("{tool}/{command}/{sub_command}/{batch}/")
    shell:
        """
        hicQC --logfiles {input}\
              --labels {params.labels}\
              --outputFolder {output}
        """

rule hicCorrelate_per_batch:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
        labels = h5PerBatchLabels,
        additional = "--plotNumbers --plotFileFormat pdf"
    threads:
        32
    input:
        files = h5PerBatchFiles
    output:
        heatmap = "hicexplorer/hicCorrelate/perBatch/{command}/{subcommand}/{batch}_heatmap.pdf",
        scatterplot = "hicexplorer/hicCorrelate/perBatch/{command}/{subcommand}/{batch}_scatterplot.pdf",
    shell:
        """
        hicCorrelate --matrices {input.files}\
                     --labels {params.labels}\
                     --threads {threads}\
                     --outFileNameHeatmap {output.heatmap}\
                     --outFileNameScatter {output.scatterplot}
                        
        """

rule hicCorrelate_per_sample:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
        labels = h5PerSampleLabels,
        additional = "--plotNumbers --plotFileFormat pdf"
    threads:
        32
    input:
        files = h5PerSampleFiles
    output:
        heatmap = "hicexplorer/hicCorrelate/perSample/{command}/{subcommand}/{sample}_heatmap.pdf",
        scatterplot = "hicexplorer/hicCorrelate/perSample/{command}/{subcommand}/{sample}_scatterplot.pdf"
    shell:
        """
        hicCorrelate --matrices {input.files}\
                     --labels {params.labels}\
                     --threads {threads}\
                     --outFileNameHeatmap {output.heatmap}\
                     --outFileNameScatter {output.scatterplot}
                        
        """

rule hicSumMatrices:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        1
    params:
    threads:
        1
    input:
        files = h5PerSampleFiles
    output:
        matrix = "hicexplorer/hicSumMatrices/{command}/{subcommand}/{sample}.h5"
    shell:
        """
            hicSumMatrices --matrices {input.files} {output.matrix}
        """
