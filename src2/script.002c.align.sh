#!/bin/bash

#SBATCH --job-name=002d
#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.002d.out.txt   
#SBATCH --error=messages/messages.002d.err.txt

#
# info
# not sure, remove this comment. it seems that in this run, one sample takes 665 min (11 hours) of elapsed at 32 threads
# pick first or last sample and run at 64, 32 and 16 threads to see the elapsed time.
# fix that messages003 is actually 002c
# 002b at 64 threads took around 769 min (13 hours) of elapsed time for the alignment. 
# save the runners xx2b, c, d, etc

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
# 2. copy inputs to scratch directory. this step takes xx
#
echo "about to copy files into scratch"
date
cp /hpcdata/Mimir/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa $tdir/.
cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_3 $tdir/.
cp /hpcdata/Mimir/adrian/research/oskjuhlid/results/001_bam/si_ctrl_3.sorted.bam* $tdir/.
date

#
# 3. run eventalign for control samples
#
echo "about to align"
date

#cd $tdir/siCtrl_1/20240301_1458_1A_PAS65957_6f46460d
#time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_ctrl_1.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS65957_6f46460d_b2b6b79a.txt --threads 64 > $tdir/eventalign.si_ctrl_1.txt

#cd $tdir/siCtrl_2/20240301_1458_1B_PAS23131_a71fd82f
#time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_ctrl_2.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS23131_a71fd82f_acf4e742.txt --threads 64 > $tdir/eventalign.si_ctrl_2.txt

cd $tdir/siCtrl_3/20240301_1458_1C_PAS23868_39694623
time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_ctrl_3.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS23868_39694623_1daaac14.txt --threads 64 > $tdir/eventalign.si_ctrl_3.txt

#cd $tdir/siFTO_1/20240301_1458_1D_PAS67735_16d64039
#time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_fto_1.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS67735_16d64039_e566cd1c.txt --threads 64 > $tdir/eventalign.si_fto_1.txt

#cd $tdir/siFTO_2/20240301_1458_1E_PAS66599_6555e094
#time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_fto_2.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS66599_6555e094_49f980cb.txt --threads 32 > $tdir/eventalign.si_fto_2.txt

#cd $tdir/siFTO_3/20240301_1458_1F_PAS67788_1beb9119
#time nanopolish eventalign --reads fastq_pass/merged.fastq.gz --bam $tdir/si_fto_3.sorted.bam --genome $tdir/Homo_sapiens.GRCh38.cdna.all.fa --scale-events --signal-index --progress --summary sequencing_summary_PAS67788_1beb9119_384e4858.txt --threads 32 > $tdir/eventalign.si_fto_3.txt

date

#
# 4. copy back results to /hpcdata
echo "about to copy results out of scratch to my dirs"
date
cp $tdir/eventalign* /hpcdata/Mimir/adrian/research/oskjuhlid/results/tmp/.
date

#
# 5. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."
date
