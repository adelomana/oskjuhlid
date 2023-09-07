#!/bin/bash

# run as: ./sender.02.mapping.sh > messages.02.txt 2>&1

date
python --version
echo "about to start working"

#
# 1. generate bam files on cDNA on control
#
echo "about to map reads"
time minimap2 -ax map-ont \
     -t 4 \
     /Users/adrian/databases/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa \
     /Users/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq \
    | samtools sort -o /Users/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam

echo "about to index generated bam file"
time samtools index /Users/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam
time samtools quickcheck /Users/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam

#
# 2. generate bam files on cDNA on FTO
#
echo "about to map reads"
time minimap2 -ax map-ont \
     -t 4 \
     /Users/adrian/databases/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa \
     /Users/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/fastq_pass/merged.fastq \
    | samtools sort -o /Users/adrian/projects/oskjuhlid/results/bam/fto_cdna.sorted.bam

echo "about to index generated bam file"
time samtools index /Users/adrian/projects/oskjuhlid/results/bam/fto_cdna.sorted.bam
time samtools quickcheck /Users/adrian/projects/oskjuhlid/results/bam/fto_cdna.sorted.bam

echo "all done."
date
