/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_hic/Snakefile ${1}\
    --configfile /home/150/sxk150/mcf10a_hic/config.yaml\
	--use-conda\
	-d ~/data/mcf10a-hic \
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_hic/cluster.json\
        --keep-going\
	-pr ${2}\
	--config rest_enzyme=DpnII_HinfI machine=gadi project=PRJEB21971 library_typ=pe
