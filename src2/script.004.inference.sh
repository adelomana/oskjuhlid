#!/bin/bash

#SBATCH --job-name=004
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.004.out.txt   
#SBATCH --error=messages/messages.004.err.txt

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
# 2. copy files into scratch
#
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/results/003_dataprep $tdir/.
date

#
# 3. process files---inference
#
echo "about to run inference"
cd $tdir/003_dataprep

date
time m6anet inference --input_dir si_ctrl_1 --out_dir inference_control_1 --n_processes 4 --num_iterations 1000
date

date
time m6anet inference --input_dir si_ctrl_2 --out_dir inference_control_2 --n_processes 4 --num_iterations 1000
date

date
time m6anet inference --input_dir si_ctrl_3 --out_dir inference_control_3 --n_processes 4 --num_iterations 1000
date

date
time m6anet inference --input_dir si_fto_1 --out_dir inference_fto_1 --n_processes 4 --num_iterations 1000
date

date
time m6anet inference --input_dir si_fto_2 --out_dir inference_fto_2 --n_processes 4 --num_iterations 1000
date

date
time m6anet inference --input_dir si_fto_3 --out_dir inference_fto_3 --n_processes 4 --num_iterations 1000
date

#
# 4. copy files back to /hpcdata
#
echo "about to copy results out of scratch to my dirs"
date
cp -rf inference* /hpcdata/Mimir/adrian/research/oskjuhlid/results/004_inference/.
date

#
# 5. clean scratch
#
echo "about to remove scratch dir"
date
rm -rf $tdir
date
echo "all done."
