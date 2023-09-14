#!/bin/bash

#SBATCH --job-name=example
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.04.out.txt   
#SBATCH --error=messages/messages.04.err.txt
#SBATCH --partition=mimir

#
# info | this job took xx hours and xx min in real elapsed time.
#

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
/usr/bin/time -v cp /hpcdata/Mimir/adrian/projects/oskjuhlid/results/alignments/eventalign.control.txt $tdir/.

#
# 3. run dataprep
#
echo "about to align"
cd $tdir
mkdir control
/usr/bin/time -v m6anet dataprep --eventalign eventalign.control.txt \
	      --out_dir control \
	      --n_processes 8

echo "about to copy results out of scratch to my dirs"
time cp -rf control /hpcdata/Mimir/adrian/projects/oskjuhlid/results/data_prep/.
cd -
pwd

#
# 4. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
