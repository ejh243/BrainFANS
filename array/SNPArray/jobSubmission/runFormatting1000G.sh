#!/bin/bash
#SBATCH -V # export all environment variables to the batch job.
#SBATCH -q sq # submit to the serial queue
#SBATCH -l walltime=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH -e format1KG.err # error file
#SBATCH -o format1KG.log # output file


## print start date and time
echo Job started on:
date -u

source $1 || exit

cd ${SCRIPTDIR}/array/SNPArray/preprocessing/utilitys

module load BCFtools/1.9-foss-2018b

sh Download1000GData.sh
sh Format1000GPlinkgr38.sh

## print finish date and time
echo Job finished on:
date -u

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/format1KG_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/format1KG_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/format1KG_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/format1KG_${SLURM_JOB_ID}.err"