#!/bin/bash

#SBATCH --partition=mimir         # request node from a specific partition
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=2       # multiply this number by the number of nodes to know the number of threads requested
#SBATCH --hint=nomultithread      # Suppress multithread
#SBATCH --output=messages/messages.000.out.txt   
#SBATCH --error=messages/messages.000.err.txt

date
echo "about to start working"

# 1. transfer files
time rsync -azPvh /hpcdata/Mimir/shared/kak101_Nanopore2/siCtrl_1 /hpcdata/Mimir/adrian/research/oskjuhlid/data/.
time rsync -azPvh /hpcdata/Mimir/shared/kak101_Nanopore2/siCtrl_2 /hpcdata/Mimir/adrian/research/oskjuhlid/data/.
time rsync -azPvh /hpcdata/Mimir/shared/kak101_Nanopore2/siCtrl_3 /hpcdata/Mimir/adrian/research/oskjuhlid/data/.
time rsync -azPvh /hpcdata/Mimir/shared/kak101_Nanopore2/siFTO_1 /hpcdata/Mimir/adrian/research/oskjuhlid/data/.
time rsync -azPvh /hpcdata/Mimir/shared/kak101_Nanopore2/siFTO_2 /hpcdata/Mimir/adrian/research/oskjuhlid/data/.
time rsync -azPvh /hpcdata/Mimir/shared/kak101_Nanopore2/siFTO_3 /hpcdata/Mimir/adrian/research/oskjuhlid/data/.

echo "all done."
date
cat /proc/cpuinfo
