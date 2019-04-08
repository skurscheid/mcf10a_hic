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

##### load additional functions #####

include: "scripts/helper.py"

##### load additional workflow rules #####

include: "rules/fastp.smk"

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
                file = fastp_targets(units)]),
        expand("fastp/report/se/{file}.end1.fastp.{suffix}",
                file = fastp_targets(units),
                suffix = ["json", "html"]),
        expand("fastp/report/se/{file}.end2.fastp.{suffix}",
                file = fastp_targets(units),
                suffix = ["json", "html"])


rule all_align:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

