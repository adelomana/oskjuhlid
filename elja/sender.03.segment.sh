#!/bin/bash

#SBATCH --partition=mimir             # request node from a specific partition
#SBATCH --nodes=1                     # number of nodes
#SBATCH --ntasks-per-node=4          # multiply this number by the number of nodes to know the number of threads requested
#SBATCH --hint=nomultithread          # Suppress multithread
#SBATCH --output=messages.03.out.txt   
#SBATCH --error=messages.03.err.txt

date
python --version
echo "about to start working"

# Location of scratch directory on the compute nodes
scratchlocation=/scratch/users

# Create a user directory if it does not exist
if [ ! -d $scratchlocation/$USER ]; then
    mkdir -p $scratchlocation/$USER
fi

# Create a temporary directory with a unique identifier associated with your jobid
tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo $tdir

#
# 1. run eventalign from nanopolish
#
echo "about to segment"
time /users/home/adrian/software/nanopolish/nanopolish eventalign --reads /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq --bam /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam --genome /users/home/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --summary /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/sequencing_summary_PAG78117_b8f781c9_cd40db60.txt --threads 4 > eventalign.txt

#
# end of adrian work
#

# IMPORTANT. Delete the temporary directory and all of its content
rm -rf $tdir
echo "all done."
date
