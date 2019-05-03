from pathlib import Path
import os
import yaml
import pandas as pd

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
        t.append(row['batch'] + "/" + row['sample'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

def hicmatrixbuilder_targets(units):
    """function for creating snakemake targets for executing hicmatrixbuilder rule"""
    t = []
    for index, row in units.iterrows():
        t.append(row['batch'] + "/" + row['sample'] + "/" + row['sample'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

def hicQCInput(wildcards):
    """function for fetching QC log files per batch"""
    t = []
    for index, row in units[units.batch == wildcards["batch"]].iterrows():
        t.append(wildcards["tool"] + "/" + wildcards["command"] + "/" + wildcards["sub_command"] + "/" + row['batch'] + "/" + row['sample'] + "/" + row['sample'] + "_" + row['lane'] + "_" + str(row['replicate']) + "/qc/QC.log")
    return(t)

def hicQCLabels(wildcards):
    """function for fetching QC log files per batch"""
    t = []
    for index, row in units[units.batch == wildcards["batch"]].iterrows():
        t.append(row['sample'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)

