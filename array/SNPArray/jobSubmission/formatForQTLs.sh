#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=SNPArray/logFiles/SNPQC.o
#SBATCH --error=SNPArray/logFiles/SNPQC.e
#SBATCH --job-name=SNPQC


## print start date and time
echo Job started on:
date -u

source $1

sh preprocessing/recodeForQTLs.sh ${IMPUTEDIR}/ImputationOutput/All/hg38 ${QTLDIR}/All/hg38
sh preprocessing/recodeForQTLs.sh ${IMPUTEDIR}/ImputationOutput/EUR/hg38 ${QTLDIR}/EUR/hg38


## print end date and time
echo Job ended on:
date -u