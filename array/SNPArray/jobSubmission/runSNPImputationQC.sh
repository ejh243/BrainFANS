#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=QCImputation_%j.log
#SBATCH --error=QCImputation_%j.err
#SBATCH --job-name=QCImputation


## print start date and time
echo Job started on:
date -u


source $1 || exit

echo "Processing data located in :" ${DATADIR}


cd ${IMPUTEDIR}/ImputationOutput

module load $RVERS

cd ${SCRIPTDIR}/array/SNPArray/preprocessing

Rscript 6_summarizeImputation.r All/ ${KGG}/1000GP_Phase3_combined.legend ALL
Rscript 6_summarizeImputation.r EUR/ ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab AF

module purge
module load VCFtools
module load BCFtools

## combine imputation output separately for ALL and EUR versions
sh 7_combineImputationOutput.sh ${IMPUTEDIR}/ImputationOutput/All/hg38
sh 7_combineImputationOutput.sh ${IMPUTEDIR}/ImputationOutput/EUR/hg38

## reformat for use with verifyBamID
sh reformatForVerifyBamID.sh ${IMPUTEDIR}/ImputationOutput/All/hg38

## print end date and time
echo Job finished:
date -u


mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/QCImputation_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/QCImputation_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/QCImputation_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/QCImputation_${SLURM_JOB_ID}.err"