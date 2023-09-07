#!/bin/bash

date
#conda activate py37
python --version

echo "about to start working"

# 1. generate index
time nanopolish index \
     -d /Users/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fast5_pass/ \
     /Users/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq \
     -s /Users/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/sequencing_summary_PAG78117_b8f781c9_cd40db60.txt \
     --verbose


time nanopolish index \
     -d /Users/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/fast5_pass/ \
     /Users/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/fastq_pass/merged.fastq \
     -s /Users/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/sequencing_summary_PAG77644_b6ab55c1_7d92c353.txt \
     --verbose

echo "all done."
date
