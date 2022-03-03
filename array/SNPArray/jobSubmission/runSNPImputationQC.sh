#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=SNPArray/logFiles/QCImputation.o
#SBATCH --error=SNPArray/logFiles/QCImputation.e
#SBATCH --job-name=QCImputation



## print start date and time
echo Job started on:
date -u


####### 

## NOTE: Do not store confidenial information in this file use the config file

######

source $1


cd ${IMPUTEDIR}/ImputationOutput

module load R/3.6.3-foss-2020a

Rscript ${SCRIPTDIR}/preprocessing/5_summarizeImputation.r All/ ${KGG}/1000GP_Phase3_combined.legend ALL
Rscript ${SCRIPTDIR}/preprocessing/5_summarizeImputation.r EUR/ ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab AF

module purge
module load VCFtools
cd ${SCRIPTDIR}

## combine imputation output separately for ALL and EUR versions
sh preprocessing/6_combineImputationOutput.sh ${IMPUTEDIR}/ImputationOutput/All
sh preprocessing/6_combineImputationOutput.sh ${IMPUTEDIR}/ImputationOutput/EUR

## print end date and time
echo Job finished:
date -u
