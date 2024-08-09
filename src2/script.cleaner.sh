#!/bin/bash

#SBATCH --job-name=cleaner
#SBATCH --nodelist=compute-70
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.cleaner.out.txt   
#SBATCH --error=messages/messages.cleaner.err.txt

date
pwd
echo "about to start working"

tdir=/scratch/users/adrian/1178474-4DBi
rm -rf $tdir

echo "testing is gone"
ls $tdir

echo "all done."
date
