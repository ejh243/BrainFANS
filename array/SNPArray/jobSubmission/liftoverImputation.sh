#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=SNPArray/logFiles/LiftOverImputation-%A_%a.o
#SBATCH --error=SNPArray/logFiles/LiftOverImputation-%A_%a.e
#SBATCH --job-name=LiftOverImputation
#SBATCH --array=1-22




## print start date and time
echo Job started on:
date -u


####### 

## NOTE: Do not store confidenial information in this file use the config file

######

source $1

# load software modules needed
module purge
module load picard
module load BCFtools
## liftover SNPs to hg38

#sh ${SCRIPTDIR}/preprocessing/5_liftoverhg38.sh ${IMPUTEDIR}/ImputationOutput/All ${SLURM_ARRAY_TASK_ID}
sh ${SCRIPTDIR}/preprocessing/5_liftoverhg38.sh ${IMPUTEDIR}/ImputationOutput/EUR ${SLURM_ARRAY_TASK_ID}

