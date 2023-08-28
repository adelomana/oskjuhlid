#!/bin/bash

#SBATCH --partition=mimir         # request node from a specific partition
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=2       # multiply this number by the number of nodes to know the number of threads requested
#SBATCH --hint=nomultithread      # Suppress multithread
#SBATCH --output=messages/message.01.out.txt   
#SBATCH --error=messages/message.01.err.txt

date
python --version
echo "about to start working"

# 1. generate index
time /users/home/adrian/software/nanopolish/nanopolish index \
     -d /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fast5_pass/ \
     /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq \
     -s /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/sequencing_summary_PAG78117_b8f781c9_cd40db60.txt \
     --verbose


time /users/home/adrian/software/nanopolish/nanopolish index \
     -d /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/fast5_pass/ \
     /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/fastq_pass/merged.fastq \
     -s /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/sequencing_summary_PAG77644_b6ab55c1_7d92c353.txt \
     --verbose

echo "all done."
date
