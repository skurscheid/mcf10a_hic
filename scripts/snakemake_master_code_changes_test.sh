export target='all_hicBuildMatrix_bin_mcool'

/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_hic/Snakefile\
    -R `/home/150/sxk150/mcf10a_hic/scripts/cli_snakemake.sh ${target} --lc`\
    --configfile /home/150/sxk150/mcf10a_hic/config.yaml\
	--use-conda\
	--cluster "qsub -P {cluster.P}\
                    -l ncpus={threads} \
                    -q {cluster.queue} \
                    -l mem={cluster.mem} \
                    -l wd\
                    -l walltime={cluster.walltime}\
		            -l storage={cluster.storage}\
		            -l jobfs={cluster.jobfs}\
                    -M {cluster.M}\
                    -m {cluster.m}\
                    -e {cluster.error_out_dir} \
                    -o {cluster.std1_out_dir}" \
	--jobs 100\
	-d ~/data/mcf10a-hic\
	--rerun-incomplete \
    --local-cores 1\
	--cluster-config /home/150/sxk150/cellular_identitmcf10a_hicy_nucleome/cluster.json\
    --config rest_enzyme=DpnII_HinfI machine=gadi\
    --keep-going\
	-prn
