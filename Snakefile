# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"

samples = pd.read_csv(config["samples"], sep = "\t").set_index("sample", drop=False)
units = pd.read_csv(config["units"], sep = "\t").set_index(["sample", "batch", "lane", "replicate"], drop=False)

rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_trim:
    input:
    # for index, row in units.iterrows():
        #print(row['batch'], row['sample'], row['lane'], row['replicate'])
        expand("fastp/trimmed/{batch}/{sample}_{lane}_{replicate}.{end}.fastq.gz",
        "fastp/report/{batch}/{sample}_{lane}_{replicate}.fastp.{suffix}",

rule all_align:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

include: "rules/fastp.smk"
