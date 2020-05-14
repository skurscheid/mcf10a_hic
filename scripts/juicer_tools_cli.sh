#!/bin/bash

#java -Xms4g -Xmx128g -jar ~/software/bin/juicer/scripts/common/juicer_tools.jar pre $1 $2 /home/150/sxk150/data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/chrom_sizes.txt

#java -Xms4g -Xmx128g -jar ~/software/bin/juicer/scripts/common/juicer_tools.jar pre merged_nodups.txt.gz inter.hic /home/150/sxk150/data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/chrom_sizes.txt -f home/150/sxk150/software/bin/juicer/restriction_sites/hg38_Arima.txt -s inter.txt -g inter_hists.m

java -jar ~/software/bin/juicer/scripts/common/juicer_tools_1.21.01.jar hiccups -m 1024 -r 5000,10000 -c 1 --ignore-sparsity inter.hic test_loops --threads 12
java -jar ~/software/bin/juicer/scripts/common/juicer_tools_1.21.01.jar hiccups -m 1024 -r 5000,10000 --ignore-sparsity inter.hic all_chrom_loops --threads 12

# HICCUPS
# medium resolution parameters according to Rao, Huntley et al. 2014
# https://github.com/aidenlab/juicer/wiki/HiCCUPS

module load java
module load cuda
java -jar ~/software/bin/juicer/scripts/common/juicer_tools_1.21.01.jar hiccups -m 2048 -r 5000,10000,25000 -k KR -f 0.1,0.1,0.1 -p 4,2,1 -i 7,5,3 -t 0.02,1.5,1.75 -d 20000,20000,50000 --threads 12 inter.hic all_chrom_loops_medium_resolution