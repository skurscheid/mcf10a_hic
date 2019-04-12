# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####

min_version("5.1.2")

##### load codMCF10AshZsheets #####

#configfile: "config.yaml"

samples = pd.read_csv(config["samples"], sep = "\t").se(MCF10AshZ, drop=False)
units = pd.read_csv(config["units"], sep = "\t").set[MCF10AshZ, "batch", "lane", "replicate"], drop=False)

##### load additional functions #####

include: "scripts/helper.py"

##### load additional workflow rules #####

include: "rules/fastp.smk"
include: "rules/bowtie2.smk"

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

rule hicbuildmatrix_single_test_run:
    input:
        "hicexplorer/hicBuildMatrix/test_run/HindIII/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10AshZ/MCF10AshZ_L001_2_hic_matrix.h5",
        "hicexplorer/hicBuildMatrix/test_run/HindIII/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10AshZ/MCF10AshZ_L001_2/qc"

rule align_one:
    input:
        "bowtie2/align/se/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10AshZ_L001_2.end1.bam",
        "bowtie2/report/se/NB501086_0064_DTremethick_JCSMR_HiC_shZ_TGFb/MCF10AshZ_L001_2.end1.txt"
