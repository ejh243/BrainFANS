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

## format paths in config file with project name
echo "Loading config file for project: " $1
export PROJECT=$1
POP=$2
source ./DNAm/config/config.txt 

## load R config file
## default file used unless alternative provided on command line
if [ -z "$2" ];
then
    RCONFIG=$SLURM_SUBMIT_DIR/DNAm/config/config.r
else
    RCONFIG=$2
fi

module load R/3.6.3-foss-2020a

mkdir -p ${DATADIR}/4_analysis/QTLs/Input/$POP/methylation
mkdir -p ${DATADIR}/4_analysis/QTLs/Input/$POP/genotype
mkdir -p ${DATADIR}/4_analysis/QTLs/Input/$POP/covariate

Rscript DNAm/analysis/QTLs/formatInputFiles.r ${DATADIR} ${REFDIR}/EPICArray ${SNPDIR} $POP

Rscript DNAm/analysis/QTLs/formatInputFilesByCT.r ${DATADIR} ${REFDIR}/EPICArray ${SNPDIR} $POP



