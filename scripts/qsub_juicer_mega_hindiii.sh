#!/bin/bash
#PBS -P pb97
#PBS -l walltime=96:00:00
#PBS -l wd
#PBS -q biodev
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=8
#PBS -l mem=128GB
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -l storage=scratch/kv78+gdata/kv78

about=${cli_about}

module load bwa/0.7.17
module load java/jdk-8.40
~/software/bin/juicer/scripts/common/mega.sh -g hg38 -s HindIII -D ~/software/bin/juicer

