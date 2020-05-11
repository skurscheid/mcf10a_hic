#!/bin/bash

#java -Xms4g -Xmx128g -jar ~/software/bin/juicer/scripts/common/juicer_tools.jar pre $1 $2 /home/150/sxk150/data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/chrom_sizes.txt

#java -Xms4g -Xmx128g -jar ~/software/bin/juicer/scripts/common/juicer_tools.jar pre merged_nodups.txt.gz inter.hic /home/150/sxk150/data/References/Genomes/Homo_sapiens/GRCh38_ensembl84/canonical/chrom_sizes.txt -f home/150/sxk150/software/bin/juicer/restriction_sites/hg38_Arima.txt -s inter.txt -g inter_hists.m