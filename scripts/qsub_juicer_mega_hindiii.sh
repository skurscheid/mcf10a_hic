#!/bin/bash
#PBS -P pb97
#PBS -l walltime=96:00:00
#PBS -l wd
#PBS -q biodev
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=32
#PBS -l mem=190GB
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -l storage=scratch/kv78+gdata/kv78

module load java

~/software/bin/juicer/scripts/common/mega.sh -g hg38 -s HindIII -D ~/software/bin/juicer