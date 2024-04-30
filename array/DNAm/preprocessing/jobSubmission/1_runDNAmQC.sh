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

# 1. the first input is your project folder.
# 2. the second input is the full path to config file.
# 3. the third input is a the path to the config.r file.

#-----------------------------------------------------


## print start date and time
echo Job started on:
date -u
JOBNAME="QCDNAdata"

# Forcibly moves the user to the scripts directory
cd $(dirname $0)/../../../..

## needs to be executed from the scripts folder

## format paths in config file with project name
echo "Loading config file for project: " $1
export PROJECT=$1
# the second input should be config file
source $2 || exit 1
RCONFIG=$3


## load modules
module load Pandoc
#module load R/3.6.3-foss-2020a
module load $RVERS   # load specified R version
echo $RVERS

Rscript DNAm/preprocessing/loadDataGDS.r ${DATADIR}

mkdir -p ${GDSDIR}/QCmetrics

Rscript DNAm/preprocessing/calcQCMetrics.r ${DATADIR} ${REFDIR}

Rscript -e "rmarkdown::render('DNAm/preprocessing/QC.rmd', output_file='QC.html')" --args ${DATADIR} ${RCONFIG} $USER

## mv markdown report to correct location
mv DNAm/preprocessing/QC.html ${GDSDIR}/QCmetrics

Rscript DNAm/preprocessing/clusterCellTypes.r ${DATADIR} ${REFDIR}


Rscript -e "rmarkdown::render('DNAm/preprocessing/QCwithinCellType.rmd', output_file='QCwithinCellType.html')" --args ${DATADIR} ${REFDIR} $USER

## mv markdown report to correct location
mv DNAm/preprocessing/QCwithinCellType.html ${GDSDIR}/QCmetrics

Rscript DNAm/preprocessing/normalisation.r ${DATADIR} ${REFDIR}

mkdir -p ${GDSDIR}/QCmetrics/CETYGO

Rscript DNAm/preprocessing/CETYGOdeconvolution.r ${DATADIR}

## print finish date and time
echo Job finished on:
date -u