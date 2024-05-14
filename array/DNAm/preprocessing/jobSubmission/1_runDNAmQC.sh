#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=DNAm/logFiles/%u/QCDNAdata.err # error file
#SBATCH --output=DNAm/logFiles/%u/QCDNAdata.log # output file
#SBATCH --job-name=QCDNAdata

#------------------------------------------------------

# 1. the first input is the full path to config file.

#-----------------------------------------------------


## print start date and time
echo Job started on:
date -u
JOBNAME="QCDNAdata"

source $1 || exit 1

echo "Processing data located in :" ${DATADIR}

## load modules
echo "Loading R module :" $RVERS
module load Pandoc
module load $RVERS   # load specified R version

cd ${SCRIPTSDIR}/array/DNAm/preprocessing/

Rscript loadDataGDS.r ${DATADIR}

chmod 755 ${DATADIR}/2_gds/raw.gds


mkdir -p ${GDSDIR}/QCmetrics

Rscript calcQCMetrics.r ${DATADIR} ${REFDIR}

Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args ${DATADIR} ${RCONFIG} $USER

## mv markdown report to correct location
mv QC.html ${GDSDIR}/QCmetrics/

Rscript clusterCellTypes.r ${DATADIR} ${REFDIR}

Rscript -e "rmarkdown::render('QCwithinCellType.rmd', output_file='QCwithinCellType.html')" --args ${DATADIR} ${REFDIR} $USER

## mv markdown report to correct location
mv QCwithinCellType.html ${GDSDIR}/QCmetrics

Rscript normalisation.r ${DATADIR} ${REFDIR}
chmod 755 ${DATADIR}/2_gds/rawNorm.gds

mkdir -p ${GDSDIR}/QCmetrics/CETYGO

Rscript CETYGOdeconvolution.r ${DATADIR}

## print finish date and time
echo Job finished on:
date -u