#!/bin/sh
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

## print start date and time
echo Job started on:
date -u
JOBNAME="QCDNAdata"
    
## needs to be executed from the scripts folder
echo "Changing Folder to: " $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1
source ./DNAm/config/config.txt 

module load Pandoc
module load R/3.6.3-foss-2020a

Rscript DNAm/preprocessing/1_loadDataGDS.r ${DATADIR}

mkdir -p ${GDSDIR}/QCmetrics

Rscript DNAm/preprocessing/2_calcQCMetrics.r ${DATADIR} ${REFDIR}

Rscript -e "rmarkdown::render('DNAm/preprocessing/3_QC.rmd', output_file='QC.html')" --args ${DATADIR} $2 $USER

## mv markdown report to correct location
mv DNAm/preprocessing/QC.html ${GDSDIR}/QCmetrics

Rscript DNAm/preprocessing/4_clusterCellTypes.r ${DATADIR} $2


Rscript -e "rmarkdown::render('DNAm/preprocessing/5_QCwithinCellType.rmd', output_file='QCwithinCellType.html')" --args ${DATADIR} $2 $USER

## mv markdown report to correct location
mv DNAm/preprocessing/QCwithinCellType.html ${GDSDIR}/QCmetrics

## print finish date and time
echo Job finished on:
date -u