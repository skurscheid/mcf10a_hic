#!/bin/bash
__conda_setup="$('/home/skurscheid/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/skurscheid/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/skurscheid/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/skurscheid/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
~/miniconda3/envs/snakemake-testing/bin/snakemake -s $1 $2\
	--configfile ./config.yaml \
	--use-conda\
	--jobs 80\
	-d /home/skurscheid/data/mcf10a-hic\
	--rerun-incomplete \
        --local-cores 4\
	-pr\
