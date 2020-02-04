import pandas as pd
import yaml
stream = open('config.yaml', 'r')
config = yaml.load(stream, Loader=yaml.SafeLoader)

wildcards = dict()
wildcards = {"batch" : "170306_NB501086_0102_HiC1_6_run4",
             "tool" : "hicexplorer",
             "command" : "hicBuildMatrix_bin",
             "sub_command" : "10000"}

runTable = pd.read_csv("SraRunTable.txt", sep = ",").set_index(["Run", "BioProject", "BioSample"])


