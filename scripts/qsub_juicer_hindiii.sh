#!/bin/bash
#PBS -P pb97
#PBS -l walltime=48:00:00
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

~/software/bin/juicer/scripts/juicer.sh -g hg38 -a $about -p ~/data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/chrom_sizes.txt -y ~/software/bin/juicer/restriction_sites/hg38_HindIII.txt -D ~/software/bin/juicer -s HindIII -t $PBS_NCPUS

