#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=LiftOverImputation_%A_%a.log
#SBATCH --error=LiftOverImputation_%A_%a.err
#SBATCH --job-name=LiftOverImputation
#SBATCH --array=1-22




## print start date and time
echo Job started on:
date -u

source $1 || exit

# load software modules needed
module purge
module load picard
module load BCFtools

${SCRIPTDIR}/preprocessing/

## liftover SNPs to hg38

sh 5_liftoverhg38.sh ${IMPUTEDIR}/ImputationOutput/All ${SLURM_ARRAY_TASK_ID}
sh 5_liftoverhg38.sh ${IMPUTEDIR}/ImputationOutput/EUR ${SLURM_ARRAY_TASK_ID}


## print finish date and time
echo Job finished on:
date -u


mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/LiftOverImputation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${DATADIR}/logFiles/LiftOverImputation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
mv "${SLURM_SUBMIT_DIR}/LiftOverImputation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${DATADIR}/logFiles/LiftOverImputation_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"