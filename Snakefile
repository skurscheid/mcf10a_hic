# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####

min_version("5.1.2")

##### set variables #####

machine = config['machine']
runTable_file = config['params']['general'][config['project']]['runTable']['file']
uibrary_type = config['library_type']
selected_columns = config['params']['general'][config['project']]['runTable']['selected_columns']
chip_input_value = config['params']['general'][config['project']]['runTable']['chip_input_value']
rest_enzyme = config['rest_enzyme']

runTable = pd.read_csv(config['params']['general'][project]['runTable']['file'], sep = ',', index_col='row_id')


##### load additional functions #####

include: "scripts/helper.py"

##### global variables/constraints #####
wildcard_constraints:
    run="[^_]*",
    resolution="\d+"

##### build targets #####
rule all:
    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_sra_download:
    input:
        expand("raw/{file}{suffix}.fastq.gz",
               file = make_targets_from_runTable_new(runTable, library_type, selected_columns, chip_input_value),
               suffix = ['_1', '_2'])

rule all_trim:
    input:
        expand("fastp/trimmed/se/{file}{end}.fastq.gz",
                file = make_targets_from_runTable(runTable),
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]])

rule all_align_global:
    input:
        expand("bowtie2/align_global/se/{file}{end}.bam",
                file = make_targets_from_runTable(runTable),
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]]),
        expand("bowtie2/align_global/se/{file}{end}.unmap.fastq",
                file = make_targets_from_runTable(runTable),
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]])

rule all_align_local:
    input:
        expand("bowtie2/align_local/se/{file}{end}.bam",
                file = make_targets_from_runTable(runTable),
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]]),
        expand("bowtie2/align_local/se/{file}{end}.unmap.fastq",
                file = make_targets_from_runTable(runTable)[3],
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]])

rule all_merge_local_global:
    input:
        expand("samtools/merge/pe/{file}{end}.bam",
                file = make_targets_from_runTable(runTable),
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]])

rule all_combine_bam_files:
    input:
        expand("mergeSam/combine/pe/{file}.bam",
               file = make_targets_from_runTable(runTable))

rule all_hicBuildMatrix_bin_test_run:
    input:
        expand("hicexplorer/hicBuildMatrix_bin/test_run/{resolution}/{file}_hic_matrix.{ext}",
               resolution = "20000",
               file = make_targets_from_runTable(runTable),
               ext = ["h5", "bam"])

rule all_hicBuildMatrix_rest_test_run:
    input:
        expand("hicexplorer/hicBuildMatrix_rest/test_run/{res_enzyme}/{file}_hic_matrix.{ext}",
               res_enzyme = "HindIII",
               file = make_targets_from_runTable(runTable),
               ext = ["h5", "bam"])
               
rule all_hicBuildMatrix_bin:
    input:
        expand("hicexplorer/hicBuildMatrix_bin/{resolution}/{file}_hic_matrix.{ext}",
               resolution = "20000",
               file = make_targets_from_runTable(runTable),
               ext = ["h5"]),

rule all_hicBuildMatrix_bin_mcool:
    input:
        expand('hicexplorer/hicBuildMatrix_bin/multi_resolution/{file}_hic_matrix.mcool',
               file = make_targets_from_runTable(runTable))
 

rule all_hicBuildMatrix_rest:
    input:
        expand("hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{file}_hic_matrix.{ext}",
               res_enzyme = "HindIII",
               file = make_targets_from_runTable(runTable),
               ext = ["h5", "bam"])

rule all_hicQC_rest:
    input:
        expand('hicexplorer/hicQC/hicBuildMatrix_rest/HindIII/{biosample}/',
                biosample = list(pd.unique(runTable.BioSample)))

rule all_hicQC_bin:
    input:
        expand('hicexplorer/hicQC/hicBuildMatrix_bin/20000/{biosample}/',
                biosample = list(pd.unique(runTable.BioSample)))

rule all_hicCorrelate_rest:
    input:
        expand('hicexplorer/hicCorrelate/hicBuildMatrix_rest/HindIII/{biosample}_scatter.pdf',
                biosample = list(pd.unique(runTable.BioSample)))

rule all_hicCorrelate_bin:
    input:
        expand('hicexplorer/hicCorrelate/hicBuildMatrix_bin/20000/{biosample}_scatter.pdf',
                biosample = list(pd.unique(runTable.BioSample)))



##### rules for extended trial runs #####
trial_samples = ['SAMN08446098/rep1/SRR6657510', 'SAMN08446098/rep1/SRR6657511',
                 'SAMN08446098/rep1/SRR6657512', 'SAMN08446098/rep1/SRR6657513',
                 'SAMN08446098/rep1/SRR6657514', 'SAMN08446097/rep2/SRR6657515',
                 'SAMN08446097/rep2/SRR6657516', 'SAMN08446097/rep2/SRR6657517',
                 'SAMN08446097/rep2/SRR6657518', 'SAMN08446097/rep2/SRR6657519']

rule trial_hicBuildMatrix_rest:
    input:
        expand("hicexplorer/hicBuildMatrix_rest/{res_enzyme}/{file}_hic_matrix.{ext}",
               res_enzyme = "HindIII",
               file = trial_samples,
               ext = ["h5", "bam"])

rule trial_combine_bam_files:
    input:
        expand("mergeSam/combine/pe/{file}.bam",
               file = make_targets_from_runTable(runTable)[0])

rule trial_hicBuildMatrix_bin_mcool:
    input:
        expand('hicexplorer/hicBuildMatrix_bin/multi_resolution/{file}_hic_matrix.mcool',
               file = make_targets_from_runTable(runTable)[0])

##### load additional workflow rules #####
include: "rules/sra_download.smk"
include: "rules/fastp.smk"
include: "rules/align.smk"
include: "rules/hicexplorer.smk"
