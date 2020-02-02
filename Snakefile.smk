# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####

min_version("5.1.2")

##### load codMCF10AshZsheets #####

#configfile: "config.yaml"

samples = pd.read_csv(config["samples"], sep = "\t").set_index("sample_id", drop=False)
units = pd.read_csv(config["units"], sep = "\t").set_index(["sample_id", "batch", "lane", "replicate"], drop=False)

##### load additional functions #####

include: "scripts/helper.py"

##### build targets #####

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_trim:
    input:
        expand("fastp/trimmed/se/{file}.end1.fastq.gz",
                file = fastp_targets(units)),
        expand("fastp/trimmed/se/{file}.end2.fastq.gz",
                file = fastp_targets(units)),
        expand("fastp/report/se/{file}.end1.fastp.{suffix}",
                file = fastp_targets(units),
                suffix = ["json", "html"]),
        expand("fastp/report/se/{file}.end2.fastp.{suffix}",
                file = fastp_targets(units),
                suffix = ["json", "html"])

rule all_align:
    input:
        expand("bowtie2/align/se/{file}.{end}.bam",
                file = fastp_targets(units),
                end = ["end1", "end2"]),
        expand("bowtie2/report/se/{file}.{end}.txt",
                file = fastp_targets(units),
                end = ["end1", "end2"])

rule test_align_rerun:
    input:
        expand("bowtie2_rerun/align/se/{file}.{end}.bam",
                file = "NB501086_0079_DTremethick_JCSMR_HiC_run1/MCF10ACA1_L001_1",
                end = ["end1", "end2"]),
        expand("bowtie2_rerun/report/se/{file}.{end}.txt",
                file = "NB501086_0079_DTremethick_JCSMR_HiC_run1/MCF10ACA1_L001_1",
                end = ["end1", "end2"])

rule all_hicbuildmatrix_HindIII:
    input:	
        expand("hicexplorer/hicBuildMatrix/HindIII/{file}_hic_matrix.h5",
               file = hicmatrixbuilder_targets(units)),
        expand("hicexplorer/hicBuildMatrix/HindIII/{file}.bam",
               file = hicmatrixbuilder_targets(units)),
	expand("hicexplorer/hicBuildMatrix/HindIII/{file}/qc",
               file = hicmatrixbuilder_targets(units)),

rule all_hicbuildmatrix_bin:
    input:
        expand("hicexplorer/hicBuildMatrix_bin/{bin_size}/{file}_hic_matrix.h5",
               file = hicmatrixbuilder_targets(units),
               bin_size = [10000, 20000, 50000, 100000]),
        expand("hicexplorer/hicBuildMatrix_bin/{bin_size}/{file}/qc",
               file = hicmatrixbuilder_targets(units),
               bin_size = [10000, 20000, 50000, 100000])

rule all_hicSumMatrices_bin:
    input:
        expand("hicexplorer/hicSumMatrices/hicBuildMatrix_bin/{bin_size}/{sample}.h5",
               sample = samples['sample_id'].unique().tolist(),
               bin_size = [10000]),

rule all_hicSumMatrices_HindIII:
    input:
        expand("hicexplorer/hicSumMatrices/hicBuildMatrix/{bin_size}/{sample}.h5",
               sample = samples['sample_id'].unique().tolist(),
               bin_size = "HindIII"),

rule test_run_hicbuildmatrix_HindIII:
    input:
        expand("hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/test_{sample}_{lane}_{replicate}_hic_matrix.{suffix}",
               res_enzyme = "HindIII",
               batch = "NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb",
               sample = "MCF10ATGFb",
               lane = "L001",
               replicate = ["1", "2"],
               suffix = ["h5", "bam"])
    #    expand("hicexplorer/hicBuildMatrix/test_run/{res_enzyme}/{batch}/{sample}/{sample}_{lane}_{replicate}/qc",
    #            res_enzyme = "HindIII",
    #            batch = "NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb",
    #            sample = "MCF10ATGFb",
    #            lane = "L001",
    #            replicate = ["1", "2"])

rule test_run_hicbuildmatrix_bin:
    input:
        expand("hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{batch}/{sample}/test_{sample}_{lane}_{replicate}_hic_matrix.{suffix}",
               resolution = "10000",
               batch = "NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb",
               sample = "MCF10ATGFb",
               lane = "L001",
               replicate = ["1", "2"],
               suffix = ["h5", "bam"])

rule test_hicCorrelate_perSample:
    input:
        expand("hicexplorer/hicCorrelate/perSample/hicBuildMatrix/{subcommand}/{sample}_{plot}.pdf",
                subcommand = "HindIII",
                sample = "MCF10A",
                plot = ["heatmap", "scatterplot"])

rule test_hicCorrelate_perBatch:
    input:
        expand("hicexplorer/hicCorrelate/perBatch/hicBuildMatrix_bin/{subcommand}/{batch}_{plot}.pdf",
                subcommand = "100000",
                batch = "NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb",
                plot = ["heatmap", "scatterplot"])

rule all_hicQC:
    input:
        expand("hicexplorer/hicQC/HindIII/{batch}/",
               batch = ["170306_NB501086_0102_HiC1_6_run4", "NB501086_0088_DTremethick_JCSMR_HiC_shZ_TGFb", 
                        "NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb", "NB501086_0100_DTremethick_HiC1_6_run3",
                        "NB501086_0079_DTremethick_JCSMR_HiC_run1", "NB501086_0103_DTremethick_JCSMR_HiC_shZ_TGFb_run3",
                        "NB501086_0080_DTremethick_JCSMR_HiC_run2", "Project_SN877_0303_Max_Nekrasov_Human_Breast_HiC"]),
        expand("hicexplorer/hicQC/bin/{bin_size}/{batch}/",
               bin_size = [20000, 50000, 100000],
               batch = ["170306_NB501086_0102_HiC1_6_run4", "NB501086_0088_DTremethick_JCSMR_HiC_shZ_TGFb",
                        "NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb", "NB501086_0100_DTremethick_HiC1_6_run3",
                        "NB501086_0079_DTremethick_JCSMR_HiC_run1", "NB501086_0103_DTremethick_JCSMR_HiC_shZ_TGFb_run3",
                        "NB501086_0080_DTremethick_JCSMR_HiC_run2", "Project_SN877_0303_Max_Nekrasov_Human_Breast_HiC"])

##### load additional workflow rules #####
include: "rules/fastp.smk"
include: "rules/bowtie2.smk"
include: "rules/hicexplorer.smk"
