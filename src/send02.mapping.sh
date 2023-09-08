#!/bin/bash

#SBATCH --partition=mimir             # request node from a specific partition
#SBATCH --nodes=1                     # number of nodes
#SBATCH --ntasks-per-node=16          # multiply this number by the number of nodes to know the number of threads requested
#SBATCH --hint=nomultithread          # Suppress multithread
#SBATCH --output=messages/messages.02.out.txt   
#SBATCH --error=messages/messages.02.err.txt

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
# 1. generate bam files on cDNA on control
#
echo "about to map reads"
time minimap2 -ax map-ont \
     -t 16 \
     /users/home/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa \
     /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq \
    | samtools sort -o /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam \
	       -T $tdir/reads.tmp

echo "about to index generated bam file"
time samtools index /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam
time samtools quickcheck /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam

#
# 2. generate bam files on cDNA on FTO
#
echo "about to map reads"                                                                                                                                  
time minimap2 -ax map-ont \
     -t 16 \
     /users/home/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa \
     /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/fto/20230817_1451_1G_PAG77644_b6ab55c1/fastq_pass/merged.fastq \
    | samtools sort -o /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/fto_cdna.sorted.bam \
               -T $tdir/reads.tmp

echo "about to index generated bam file"
time samtools index /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/fto_cdna.sorted.bam
time samtools quickcheck /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/fto_cdna.sorted.bam

#
# end of adrian work
#

# IMPORTANT. Delete the temporary directory and all of its content
rm -rf $tdir
echo "all done."
date
