#!/bin/bash

#SBATCH --partition=mimir
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16         
#SBATCH --hint=multithread        
#SBATCH --output=messages/messages.001.out.txt   
#SBATCH --error=messages/messages.001.err.txt

# info | this script takes three hours of elapsed time

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
# 2. copy inputs to scratch directory
#
echo "copy working files to scratch"
time cp -rf /hpcdata/Mimir/adrian/database/ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa $tdir/.
time cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_1 $tdir/.
time cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_2 $tdir/.
time cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_3 $tdir/.
time cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siFTO_1 $tdir/.
time cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siFTO_2 $tdir/.
time cp -rf /hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siFTO_3 $tdir/.

#
# 3. map reads
#
echo "about to map reads"
cd $tdir
pwd
ls

time minimap2 -ax map-ont -t 16 Homo_sapiens.GRCh38.cdna.all.fa siCtrl_1/20240301_1458_1A_PAS65957_6f46460d/fastq_pass/merged.fastq.gz | samtools sort -o si_ctrl_1.sorted.bam -T $tdir/reads.tmp
time minimap2 -ax map-ont -t 16 Homo_sapiens.GRCh38.cdna.all.fa siCtrl_2/20240301_1458_1B_PAS23131_a71fd82f/fastq_pass/merged.fastq.gz | samtools sort -o si_ctrl_2.sorted.bam -T $tdir/reads.tmp
time minimap2 -ax map-ont -t 16 Homo_sapiens.GRCh38.cdna.all.fa siCtrl_3/20240301_1458_1C_PAS23868_39694623/fastq_pass/merged.fastq.gz | samtools sort -o si_ctrl_3.sorted.bam -T $tdir/reads.tmp

time minimap2 -ax map-ont -t 16 Homo_sapiens.GRCh38.cdna.all.fa siFTO_1/20240301_1458_1D_PAS67735_16d64039/fastq_pass/merged.fastq.gz | samtools sort -o si_fto_1.sorted.bam -T $tdir/reads.tmp
time minimap2 -ax map-ont -t 16 Homo_sapiens.GRCh38.cdna.all.fa siFTO_2/20240301_1458_1E_PAS66599_6555e094/fastq_pass/merged.fastq.gz | samtools sort -o si_fto_2.sorted.bam -T $tdir/reads.tmp
time minimap2 -ax map-ont -t 16 Homo_sapiens.GRCh38.cdna.all.fa siFTO_3/20240301_1458_1F_PAS67788_1beb9119/fastq_pass/merged.fastq.gz | samtools sort -o si_fto_3.sorted.bam -T $tdir/reads.tmp


echo "about to index generated bam file"
time samtools index si_ctrl_1.sorted.bam
time samtools index si_ctrl_2.sorted.bam
time samtools index si_ctrl_3.sorted.bam
time samtools index si_fto_1.sorted.bam
time samtools index si_fto_2.sorted.bam
time samtools index si_fto_3.sorted.bam

time samtools quickcheck si_ctrl_1.sorted.bam
time samtools quickcheck si_ctrl_2.sorted.bam
time samtools quickcheck si_ctrl_3.sorted.bam
time samtools quickcheck si_fto_1.sorted.bam
time samtools quickcheck si_fto_2.sorted.bam
time samtools quickcheck si_fto_3.sorted.bam

#
# 4. copy results back to /hpcdata
#
echo "about to copy results out of scratch to my dirs"
time cp -rf *.sorted.bam* /hpcdata/Mimir/adrian/research/oskjuhlid/results/001_bam/.
cd -
pwd

#
# 5. clean scratch
#
echo "about to remove scratch dir"
rm -rf $tdir
echo "all done."

date
