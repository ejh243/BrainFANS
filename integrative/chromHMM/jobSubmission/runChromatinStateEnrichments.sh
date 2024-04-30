#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/overlapChromatinAnnotations5hmc.err # error file
#SBATCH --output=LogFiles/overlapChromatinAnnotations5hmc.log # output file
#SBATCH --job-name=overlapChromatinAnnotations

## needs to be executed from the scripts folder


## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######




## generate QC metrics
module purge
module load R/3.6.3-foss-2020a
#Rscript ATACSeq/overlapChromatinAnnotations.r
Rscript hydroxy/CGEX/overlapChromatinAnnotations.r
