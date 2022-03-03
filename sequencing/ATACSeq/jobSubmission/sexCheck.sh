#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/%u/sexCheck-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/%u/sexCheck-%A_%a.e
#SBATCH --job-name=sexCheck



## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ATACSeq/config/config.txt 
echo "Project directory is: " $DATADIR


## check script directory
echo 'Script directory is: ' ${SCRIPTDIR}

## call peaks for sex chromosomes & do read counts in these peaks

module load MACS2
module load BEDTools



sh ./ATACSeq/preprocessing/12_sexChrPeaks.sh

module load R/3.6.3-foss-2020a
Rscript ATACSeq/preprocessing/13_collateSexChecks.r ${DATADIR}/

