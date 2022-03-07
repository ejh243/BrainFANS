#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=4 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/EWAS_Sims_1.err # error file
#SBATCH --output=LogFiles/EWAS_Sims_1.log # output file
#SBATCH --job-name=EWASSimulations

## Output some useful job information

echo SBATCH: working directory is $SLURM_SUBMIT_DIR 
echo SBATCH: job identifier is $SLURM_JOBID

## print start date and time
echo Job started on:
date -u

cd $SLURM_SUBMIT_DIR

module load R/3.6.0-foss-2019a

Rscript DNAm/EWAS/Simulations/NullSimulationsAll.r DNAm/rmdConfig.mrc 5

Rscript DNAm/EWAS/Simulations/NullSimulationsAll.r DNAm/rmdConfig.mrc 6
