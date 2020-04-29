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

rule hicBuildMatrix_bin:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        4
    params:
        inputBufferSize = 400000
    threads:
        16
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

rule hicBuildMatrix_bin_mcool:
    conda:
        "../envs/hicexplorer.yaml"
    version:
        4
    params:
        inputBufferSize = 800000,
        resolution = ['20000', '40000', '100000', '1000000', '2000000'],
        genomeAssembly = config['params']['hicexplorer']['hicBuildMatrix']['genomeAssembly'],
        danglingSequence = config['params']['hicexplorer']['hicBuildMatrix']['danglingSequence']
    threads:
        16
    input:
        mate1 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end1_suffix"] + ".bam",
        mate2 = lambda wildcards: "/".join(["samtools", "sort", "se", wildcards["biosample"], wildcards["rep"], wildcards["run"]]) + config["params"]["general"]["end2_suffix"] + ".bam"
    log:
        "logs/hicexplorer/hicBuildMatrix_bin/multi_resolution/{biosample}/{rep}/{run}_mcool.txt"
    benchmark:
        "benchmarks/hicexplorer/hicBuildMatrix_bin/multi_resolution/{biosample}/{rep}/{run}_mcool.tsv"
    output:
        outHicMatrix = "hicexplorer/hicBuildMatrix_bin/multi_resolution/{biosample}/{rep}/{run}_hic_matrix.mcool"
    shell:
        """
            hicBuildMatrix --samFiles {input.mate1} {input.mate2} \
                    --threads {threads} \
                    --binSize {params.resolution}\
                    --inputBufferSize {params.inputBufferSize} \
                    --genomeAssembly {params.genomeAssembly} \
                    --danglingSequence {params.danglingSequence} \
                    --outFileName {output.outHicMatrix} 1>{log} 2>{log}
        """


rule hicQC:
    conda:
        '../envs/hicexplorer.yaml'
    version:
        1
    params:
        inputBufferSize = 400000,
        labels = lambda wildcards: list(runTable.loc[runTable.BioSample == wildcards['biosample']].Run)
    threads:
        16
    log:
        logfile = 'logs/hicexplorer/hicQC/{tool}/{res}/{biosample}.log'
    input:
        qc_files = lambda wildcards: expand('/'.join(['hicexplorer', wildcards['tool'], wildcards['res'], wildcards['biosample']]) +\
                                            '/{run}/qc/QC.log',\
                                            run = list(runTable.loc[runTable.BioSample == wildcards['biosample'], ['replicate', 'Run']].apply(lambda x: '/'.join(x), axis=1)))
    output:
        directory('hicexplorer/hicQC/{tool}/{res}/{biosample}/')
    shell:
        '''
            hicQC --logfiles {input.qc_files}\
                  --labels {params.labels}\
                  --outputFolder {output} 2>{log.logfile}
        '''

rule hicCorrelate:
    conda:
        '../envs/hicexplorer.yaml'
    version:
        1
    params:
        inputBufferSize = 400000,
        labels = lambda wildcards: list(runTable.loc[runTable.BioSample == wildcards['biosample']].Run),
        method = 'pearson',
        plotFileFormat = 'pdf',
        cli_params = '--log1p ',
        range = '20000:2000000'
    threads:
        16
    log:
        logfile = 'logs/hicexplorer/hicCorrelate/{tool}/{res}/{biosample}.log'
    input:
        matrix = lambda wildcards: expand('/'.join(['hicexplorer', wildcards['tool'], wildcards['res'], wildcards['biosample']]) +\
                                            '/{run}_hic_matrix.h5',\
                                            run = list(runTable.loc[runTable.BioSample == wildcards['biosample'], ['replicate', 'Run']].apply(lambda x: '/'.join(x), axis=1)))
    output:
        scatter = 'hicexplorer/hicCorrelate/{tool}/{res}/{biosample}_scatter.pdf',
        heatmap = 'hicexplorer/hicCorrelate/{tool}/{res}/{biosample}_heatmap.pdf'
    shell:
        '''
            hicCorrelate --matrices {input.matrix}\
                         --method {params.method}\
                         --labels {params.labels}\
                         --outFileNameHeatmap {output.heatmap}\
                         --outFileNameScatter {output.scatter}\
                         {params.cli_params}\
                         --range {params.range}\
                         --threads {threads} 2>{log.logfile}
        '''
