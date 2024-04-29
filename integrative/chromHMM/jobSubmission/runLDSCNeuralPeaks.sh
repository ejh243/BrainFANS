#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/LDSCNeuralPeaks.err # error file
#SBATCH --output=LogFiles/LDSCNeuralPeaks.log # output file
#SBATCH --job-name=LDSCNeuralPeaks

## needs to be executed from the scripts folder


## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######



module purge
module load BEDOPS
module load BEDTools
module load R/3.6.3-foss-2020a
#Rscript /gpfs/mrc0/projects/Research_Project-MRC190311/scripts/LDScoreRegression/createAnnotationFiles.r


module purge
module load Anaconda3/2020.02
source activate ldsc

source LDScoreRegression/config.txt
source LDScoreRegression/runPartionedHeritabilityOnPeaks.sh

## print start date and time
echo Job finished at:
date -u