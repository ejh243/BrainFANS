#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=240:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/bedtobam-%A_%a.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/bedtobam-%A_%a.e
#SBATCH --job-name=bedtobam

## print start date and time
echo Job started on:
date -u

start=`date +%s`

## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT

PROJECT=$2

source ./integrative/chromHMM/config/config.txt

cd $DATADIR/1_raw
FILES=($(ls *.gz))

DATA=$(echo ${FILES[${SLURM_ARRAY_TASK_ID}]%.tagAlign.gz} )
DATA=$(echo ${DATA#E073-} )

source ${SCRIPTDIR}/integrative/chromHMM/config/config.txt

mkdir -p $INDIR
echo $COMPDIR
echo
#-----------------------------------------------------------------------#

module load BEDTools
module load SAMtools

zcat ${FILES[${SLURM_ARRAY_TASK_ID}]} | bedToBam -i - -g /lustre/projects/Research_Project-MRC190311/references/hg19/hg19.chrom.sizes > $INDIR/downsampled.0.bam


## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}