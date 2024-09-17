#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=CTEWASSIMS-%A_%a.err # error file
#SBATCH --output=CTEWASSIMS-%A_%a.log # output file
#SBATCH --job-name=CTEWASSIMS
#SBATCH --array=1-10

#------------------------------------------------------

# 1. the first input is the full path to config file.

#-----------------------------------------------------

## print start date and time
echo Job started on:
date -u

echo Log files intially stored in: ${SLURM_SUBMIT_DIR}/CTEWASSIMS_${SLURM_JOB_ID}.log and ${SLURM_SUBMIT_DIR}/CTEWASSIMS_${SLURM_JOB_ID}.err

source $1 || exit 1

echo "Processing data located in :" ${DATADIR}

echo "Random seed: " ${SLURM_ARRAY_TASK_ID}

module load R/3.6.3-foss-2020a

cd ${SCRIPTSDIR}/array/DNAm/analysis/methodsDevelopment/

meanDiff=0.05
Rscript simulateCellSpecificEWAS.r ${DATADIR} ${SLURM_ARRAY_TASK_ID} ${meanDiff}

# print end date and time
echo Job finished on:
date -u
   

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/CTEWASSIMS_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log" \
"${DATADIR}/logFiles/CTEWASSIMS_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
mv "${SLURM_SUBMIT_DIR}/CTEWASSIMS_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" \
"${DATADIR}/logFiles/CTEWASSIMS_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"

echo Log files moved to: ${DATADIR}/logFiles/   
