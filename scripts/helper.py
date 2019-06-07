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

class hicCorrelateParams:
    """ functions for providing inputs and parameters for hicCorrelate """
    def __init__(self):
        self.data = []
        self.batch = []
        self.sample = []
        self.labels = []


    def h5PerBatch(self, units, wildcards):
        """function for fetching h5 matrix files per batch"""
        for index, row in units[units.batch == wildcards["batch"]].iterrows():
            self.labels.append(row['sample_id'])
            self.batch.append("/".join(["hicexplorer/hicBuildMatrix",
                                        wildcards["subcommand"],
                                        row['batch'],
                                        row['sample_id'],
                                        row['sample_id']]) + "_" + "_".join([row['lane'], str(row['replicate']), "hic_matrix.h5"]))
        
    def h5PerSample(self, units, wildcards):
        """function for fetching h5 matrix files per sample"""
        for index, row in units[units.sample_id == wildcards["sample"]].iterrows():
            self.labels.append(row['batch'])
            self.sample.append("/".join(["hicexplorer/hicBuildMatrix",
                                         wildcards["subcommand"],
                                         row['batch'],
                                         row['sample_id'],
                                         row['sample_id']]) + row['lane'] + "_" + str(row['replicate']) + "_hic_matrix.h5")

