# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.

import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####

min_version("5.1.2")

##### load codMCF10AshZsheets #####

#configfile: "config.yaml"

runTable = pd.read_csv("SraRunTable.txt", sep = ",")

##### load additional functions #####

include: "scripts/helper.py"

##### build targets #####

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_trim:
    input:
        expand("fastp/trimmed/se/{file}_{end}.fastq.gz",
                file = list(runTable["Run"])[0],
                end = ["1", "2"]),
        expand("fastp/report/se/{file}_{end}.fastp.{suffix}",
                file = list(runTable["Run"])[0],
                end = ["1", "2"],
                suffix = ["json", "html"])

rule all_align:
    input:
        expand("bowtie2/align/se/{file}_{end}.bam",
                file = list(runTable["Run"])[0],
                end = ["1", "2"]),
        expand("bowtie2/report/se/{file}_{end}.txt",
                file = list(runTable["Run"])[0],
                end = ["1", "2"]),

##### load additional workflow rules #####
include: "rules/fastp.smk"
include: "rules/bowtie2.smk"
include: "rules/hicexplorer.smk"
