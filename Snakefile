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


rule test_hicbuildmatrix_HindIII:
    input:
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_1_hic_matrix.h5",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_1.bam",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_2_hic_matrix.h5",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_2.bam",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_1/qc",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_2/qc"



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

rule hicbuildmatrix_100k_resolution_test:
    input:
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_1_hic_matrix.h5",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_2_hic_matrix.h5",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_1/qc",
        "hicexplorer/hicBuildMatrix_bin/100000/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10ATGFb/MCF10ATGFb_L001_2/qc"

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
