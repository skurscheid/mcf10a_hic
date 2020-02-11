export cli_target=${1}
/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/cellular_identity_nucleome/Snakefile ${1}\
    --configfile /home/150/sxk150/cellular_identity_nucleome/config.yaml\
	--use-conda\
	-d /scratch/kv78/cellular_identity_nucleome\
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/cellular_identity_nucleome/cluster.json\
        --keep-going\
	-pr\
	-R `/home/150/sxk150/cellular_identity_nucleome/scripts/cli_snakemake.sh all_combine_bam_files --lc`\
	${2}
