def create_testing_input(units):
    """creates test files for snakemake run"""
    for index, row in units.iterrows():
        fq1, fq2 = Path("/Users/u1001407/data/breast_cancer_testing/" ,row['fq1']), Path("/Users/u1001407/data/breast_cancer_testing/" ,row['fq2'])
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
