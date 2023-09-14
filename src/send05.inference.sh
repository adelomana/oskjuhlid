#!/bin/bash

#SBATCH --job-name=inference
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.05.out.txt   
#SBATCH --error=messages/messages.05.err.txt
#SBATCH --partition=mimir

#
# info
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
time cp -rf /hpcdata/Mimir/adrian/projects/oskjuhlid/results/data_prep/control $tdir/.
time cp -rf /hpcdata/Mimir/adrian/projects/oskjuhlid/results/data_prep/fto $tdir/.

#
# 3. run dataprep
#
echo "about to execute inference from m6anet"
cd $tdir
time m6anet inference --input_dir control --out_dir inference_control --n_processes 4 --num_iterations 1000
time m6anet inference --input_dir fto --out_dir inference_fto --n_processes 4 --num_iterations 1000

echo "about to copy results out of scratch to my dirs"
time cp -rf inference_* /hpcdata/Mimir/adrian/projects/oskjuhlid/results/inference/.
cd -
pwd

#
# 4. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
