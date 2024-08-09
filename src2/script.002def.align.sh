#!/bin/bash

#SBATCH --job-name=002def
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.002def.out.txt   
#SBATCH --error=messages/messages.002def.err.txt

date
pwd
echo "about to start working"

#
# 1. create a temporary directory with a unique identifier associated with your jobid
# 
scratchlocation=/scratch/users

if [ ! -d $scratchlocation/$USER ]; then
    mkdir -p $scratchlocation/$USER
fi

tdir=$(mktemp -d $scratchlocation/$USER/$SLURM_JOB_ID-XXXX)
echo "the scratch dir is:"
echo $tdir

#
# 2. iterate over three samples
#
cp /hpcdata/Mimir/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa $tdir/.

# 2.1. sample d
echo "sample d"
# copy inputs into scratch directory 
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siFTO_1 $tdir/.
cp /hpcdata/Mimir/adrian/research/oskjuhlid/results/001_bam/si_fto_1.sorted.bam* $tdir/.
date
# run eventalign
echo "about to align"
date
cd $tdir/siFTO_1/20240301_1458_1D_PAS67735_16d64039
time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_fto_1.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS67735_16d64039_e566cd1c.txt --threads 64 > $tdir/eventalign.si_fto_1.txt
date
# copy data out of scratch
echo "about to copy results out of scratch to my dirs"
date
cp $tdir/eventalign* /hpcdata/Mimir/adrian/research/oskjuhlid/results/tmp/.
date
# remove big files in scratch in order to continue calculating
rm -rf $tdir/si*
rm -rf $tdir/eventalign*

# 2.2. sample e
echo "sample e"
# copy inputs into scratch directory
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siFTO_2 $tdir/.
cp /hpcdata/Mimir/adrian/research/oskjuhlid/results/001_bam/si_fto_2.sorted.bam* $tdir/.
date
# run eventalign                                                                                                                                                                    
echo "about to align"
date
cd $tdir/siFTO_2/20240301_1458_1E_PAS66599_6555e094
time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_fto_2.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS66599_6555e094_49f980cb.txt --threads 32 > $tdir/eventalign.si_fto_2.txt
date
# copy data out of scratch
echo "about to copy results out of scratch to my dirs"
date
cp $tdir/eventalign* /hpcdata/Mimir/adrian/research/oskjuhlid/results/tmp/.
date
# remove big files in scratch in order to continue calculating
rm -rf $tdir/si*
rm -rf $tdir/eventalign*

# 2.3. sample f
echo "sample f"
# copy inputs into scratch directory
echo "about to copy files into scratch"
date
cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siFTO_3 $tdir/.
cp /hpcdata/Mimir/adrian/research/oskjuhlid/results/001_bam/si_fto_3.sorted.bam* $tdir/.
date
# run eventalign
echo "about to align"
date
cd $tdir/siFTO_3/20240301_1458_1F_PAS67788_1beb9119
time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_fto_3.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS67788_1beb9119_384e4858.txt --threads 32 > $tdir/eventalign.si_fto_3.txt
date
# copy data out of scratch                                                                                                                                                                                                              
echo "about to copy results out of scratch to my dirs"
date
cp $tdir/eventalign* /hpcdata/Mimir/adrian/research/oskjuhlid/results/tmp/.
date
# remove big files in scratch in order to continue calculating                                                                                                                                                                          
rm -rf $tdir/si*
rm -rf $tdir/eventalign*

#
# 3. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
