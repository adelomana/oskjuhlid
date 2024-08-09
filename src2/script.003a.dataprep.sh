#!/bin/bash

#SBATCH --job-name=003a
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.003a.out.txt   
#SBATCH --error=messages/messages.003a.err.txt

#
# info | this job took xx hours and xx min in real elapsed time.
#

date
pwd
echo "about to start working"

#
# 1. Create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users
if [ ! -d $scratchlocation/$USER ]; then
    mkdir -p $scratchlocation/$USER
fi
tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

#
# 2.1. work for sample a
#
echo "sample a"
echo "about to copy files into scratch"
date
cp /hpcdata/Mimir/adrian/research/oskjuhlid/results/002_alignments/eventalign.si_ctrl_1.txt $tdir/.
date

echo "about to run data prep"
cd $tdir
mkdir outdir
date
m6anet dataprep --eventalign eventalign.si_ctrl_1.txt --out_dir outdir --n_processes 16
date

echo "about to copy results out of scratch to my dirs"
date
cp -rf outdir /hpcdata/Mimir/adrian/research/oskjuhlid/results/003_dataprep/si_ctrl_1
date

# remove what I have in scratch
rm -rf $tdir/*

#
# 3. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date