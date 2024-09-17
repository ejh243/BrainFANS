#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=QCDNAdata_%j.err # error file
#SBATCH --output=QCDNAdata_%j.log # output file
#SBATCH --job-name=QCDNAdata

#------------------------------------------------------

# 1. the first input is the full path to config file.

#-----------------------------------------------------


## print start date and time
echo Job started on:
date -u

echo Log files intially stored in: ${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.log and ${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.err

source $1 || exit 1

echo "Processing data located in :" ${DATADIR}

## load modules
echo "Loading R module :" $RVERS
module load Pandoc
module load $RVERS   # load specified R version

cd ${SCRIPTSDIR}/array/DNAm/preprocessing/

Rscript checkRconfigFile.r ${DATADIR}

Rscript installLibraries.r ${DATADIR}

Rscript checkColnamesSampleSheet.r ${DATADIR}

mkdir -p ${GDSDIR}/QCmetrics

Rscript loadDataGDS.r ${DATADIR}

chmod 755 ${DATADIR}/2_gds/raw.gds

Rscript calcQCMetrics.r ${DATADIR} ${REFDIR}

Rscript clusterCellTypes.r ${DATADIR} ${REFDIR}

Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args ${DATADIR} ${REFDIR} ${RCONFIG} $USER

## mv markdown report to correct location
mv QC.html ${GDSDIR}/QCmetrics/

mkdir -p ${DATADIR}/3_normalised

Rscript normalisation.r ${DATADIR} ${REFDIR}
chmod 755 ${DATADIR}/2_gds/rawNorm.gds

mkdir -p ${GDSDIR}/QCmetrics/CETYGO

Rscript CETYGOdeconvolution.r ${DATADIR}

## print finish date and time
echo Job finished on:
date -u

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/QCDNAdata_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/QCDNAdata_${SLURM_JOB_ID}.err"

echo Log files moved to: ${DATADIR}/logFiles/