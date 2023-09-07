#!/bin/bash

# run as: ./sender.03.segment.sh > messages.03.txt 2>&1

# this is a very long step, computationally quite demanding. Typically the bottleneck of the whole pipeline.
# this step took xx hours @ 4 threads in an M1 chip for the control sample.
# this step generates an output of xx GB for just the control sample.

# consider running 8 threads and pray for not freezing the computer.

date
conda activate py37
python --version
echo "about to start working"

#
# 1. run eventalign from nanopolish
#
echo "about to segment control"
time nanopolish eventalign \
      --reads /Users/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq \
      --bam /Users/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam \
      --genome /Users/adrian/databases/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa \
      --scale-events --signal-index --progress \
      --summary /Users/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/sequencing_summary_PAG78117_b8f781c9_cd40db60.txt \
      --threads 6 > /Users/adrian/projects/oskjuhlid/results/alignments/eventalign.control.txt
