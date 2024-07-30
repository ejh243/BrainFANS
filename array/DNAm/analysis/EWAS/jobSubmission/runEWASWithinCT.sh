#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=EWASWithinCT_%A_%a.err # error file
#SBATCH --output=EWASWithinCT_%A_%a.log # output file
#SBATCH --job-name=SCZEWASWithinCT
#SBATCH --array=0-2

# print start date and time
echo Job started on:
date -u

## provide path to config file on submission
source $1 || exit 1

CELLTYPES=(NeuN+ Sox10+ Double-) 

module load R/4.2.1-foss-2022a

cd ${SCRIPTSDIR}/array/DNAm/analysis/EWAS/

Rscript lmWithinCT.r ${DATADIR} ${CELLTYPES[${SLURM_ARRAY_TASK_ID}]}


# print end date and time
echo Job finished on:
date -u


mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/EWASWithinCT_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${DATADIR}/logFiles/EWASWithinCT_${SLURM_JOB_ID}_${CELLTYPES[${SLURM_ARRAY_TASK_ID}]}.log"
mv "${SLURM_SUBMIT_DIR}/EWASWithinCT_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${DATADIR}/logFiles/EWASWithinCT_${SLURM_JOB_ID}_${CELLTYPES[${SLURM_ARRAY_TASK_ID}]}.err"
