#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=reformatForQTLs_%j.log
#SBATCH --error=reformatForQTLs_%j.err
#SBATCH --job-name=reformatForQTLs


## print start date and time
echo Job started on:
date -u

source $1|| exit

echo "Processing data located in :" ${DATADIR}

cd ${SCRIPTDIR}/array/SNPArray/preprocessing

sh recodeForQTLs.sh ${IMPUTEDIR}/ImputationOutput/All/hg38 ${QTLDIR}/All/hg38
sh recodeForQTLs.sh ${IMPUTEDIR}/ImputationOutput/EUR/hg38 ${QTLDIR}/EUR/hg38


## print end date and time
echo Job ended on:
date -u

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/reformatForQTLs_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/reformatForQTLs_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/reformatForQTLs_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/reformatForQTLs_${SLURM_JOB_ID}.err"