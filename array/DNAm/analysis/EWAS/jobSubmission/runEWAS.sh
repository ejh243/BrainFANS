#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=DNAm/logFiles/%u/EWAS.err # error file
#SBATCH --output=DNAm/logFiles/%u/EWAS.log # output file
#SBATCH --job-name=EWAS


# print start date and time
echo Job started on:
date -u
    
# needs to be executed from the scripts folder

# load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1
source ./DNAm/config/config.txt 

module load R/3.6.3-foss-2020a

Rscript DNAm/analysis/EWAS/lm.r ${DATADIR}
Rscript DNAm/analysis/EWAS/mlm.r ${DATADIR}
Rscript DNAm/analysis/EWAS/crr.r ${DATADIR}

# print end date and time
echo Job finished on:
date -u
    