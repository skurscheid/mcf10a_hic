from pathlib import Path
import os
import yaml
import pandas as pd

def make_targets_from_runTable(runTable, library_type, selected_columns, chip_input_value):
    t = []
    for index, row in runTable.iterrows():
        library = row[selected_columns[0]].split()[0]
        if library == chip_input_value:
            library = 'Input'
        e = list([row[selected_columns[2]], library, library_type, row['Run']])
        p = "/".join(e)
        t.append(p)
    return(t)

def create_testing_input(base_path, units):
    """creates test files for snakemake run"""
    for index, row in units.iterrows():
        fq1, fq2 = Path(base_path ,row['fq1']), Path(base_path ,row['fq2'])
        p = Path(os.path.split(fq1)[0])
        p.mkdir(parents = True, exist_ok = True)
        fq1.touch(exist_ok = True)
        fq2.touch(exist_ok = True)

def fastp_targets(units):
    """function for creating snakemake targets for executing fastp rule"""
    t = []
    for index, row in units.iterrows():
        t.append(row['batch'] + "/" + row['sample_id'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

def hicmatrixbuilder_targets(units):
    """function for creating snakemake targets for executing hicmatrixbuilder rule"""
    t = []
    for index, row in units.iterrows():
        t.append(row['batch'] + "/" + row['sample_id'] + "/" + row['sample_id'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

def hicQCInput(wildcards):
    """function for fetching QC log files per batch"""
    t = []
    for index, row in units[units.batch == wildcards["batch"]].iterrows():
        t.append(wildcards["tool"] + "/" + wildcards["command"] + "/" + wildcards["subcommand"] + "/" + row['batch'] + "/" + row['sample_id'] + "/" + row['sample_id'] + "_" + row['lane'] + "_" + str(row['replicate']) + "/qc/QC.log")
    return(t)

def hicQCLabels(wildcards):
    """function for fetching QC log files per batch"""
    t = []
    for index, row in units[units.batch == wildcards["batch"]].iterrows():
        t.append(row['sample_id'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

def h5PerBatchFiles(wildcards):
    """function for fetching h5 matrix files per batch"""
    files = []
    for index, row in units[units.batch == wildcards["batch"]].iterrows():
        files.append("/".join(["hicexplorer",
                                wildcards["command"],
                                wildcards["subcommand"],
                                row['batch'],
                                row['sample_id'],
                                row['sample_id']]) + "_" + "_".join([row['lane'], str(row['replicate']), "hic_matrix.h5"]))
    return(files)

def h5PerBatchLabels(wildcards):
    """function for fetching h5 matrix files per batch"""
    labels = []
    for index, row in units[units.batch == wildcards["batch"]].iterrows():
        labels.append(row['sample_id'])
    return(labels)

def h5PerSampleFiles(wildcards):
    """function for fetching h5 matrix files per sample"""
    files = []
    for index, row in units[units.sample_id == wildcards["sample"]].iterrows():
        files.append("/".join(["hicexplorer",
                                wildcards["command"],
                                wildcards["subcommand"],
                                row['batch'],
                                row['sample_id'],
                                row['sample_id']]) + "_" + row['lane'] + "_" + str(row['replicate']) + "_hic_matrix.h5")
    return(files)

def h5PerSampleLabels(wildcards):
    """function for fetching h5 matrix files per sample"""
    labels = []
    for index, row in units[units.sample_id == wildcards["sample"]].iterrows():
        labels.append(row['batch'])
    return(sample)

