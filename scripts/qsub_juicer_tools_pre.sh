#!/bin/bash
#PBS -P pb97
#PBS -l walltime=96:00:00
#PBS -l wd
#PBS -q biodev
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=4
#PBS -l mem=32GB
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -l storage=scratch/kv78+gdata/kv78

export about=${cli_about}
export juiceDir='/home/150/sxk150/software/bin/juicer'
export outputdir='aligned'
export genomeID='hg38'

echo $(pwd)

module load bwa/0.7.17
module load java/jdk-8.40
${juiceDir}/scripts/common/juicer_tools pre -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}

