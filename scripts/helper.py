from pathlib import Path
import os
import yaml
import pandas as pd

def make_targets_from_runTable(runTable):
    t = []
    for index, row in runTable.iterrows():
        e = list(row[['BioSample', 'replicate', 'Run']])
        p = "/".join(e)
        t.append(p)
    return(t)

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
