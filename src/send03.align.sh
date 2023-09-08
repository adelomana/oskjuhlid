#!/bin/bash

#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.03.out.txt   
#SBATCH --error=messages/messages.03.err.txt

date
python --version
pwd
echo "about to start working"

#
# 0.1. Create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users

if [ ! -d $scratchlocation/$USER ]; then
    mkdir -p $scratchlocation/$USER
fi

tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

#
# 2. copy inputs to scratch directory
#
time cp /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/fastq_pass/merged.fastq* $tdir/.
time cp /hpcdata/Mimir/adrian/projects/oskjuhlid/results/bam/control_cdna.sorted.bam* $tdir/.
time cp /hpcdata/Mimir/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa $tdir/.
time cp /hpcdata/Mimir/adrian/projects/oskjuhlid/data/august23/ctrl/20230817_1451_2G_PAG78117_b8f781c9/sequencing_summary_PAG78117_b8f781c9_cd40db60.txt $tdir/.


#
# 3. run eventalign for control samples
#
echo "about to align"
cd $tdir
time nanopolish eventalign --reads merged.fastq --bam control_cdna.sorted.bam --genome Homo_sapiens.GRCh38.cdna.all.fa \
     --scale-events --signal-index --progress \
     --summary sequencing_summary_PAG78117_b8f781c9_cd40db60.txt \
     --threads 64 > $tdir/eventalign.control.txt

echo "about to copy results out of scratch to my dirs"
time cp eventalign.control.txt /hpcdata/Mimir/adrian/projects/oskjuhlid/results/alignments/.
cd -
pwd

#
# 4. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
