#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=140:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e ../LogFiles/EWAS_LM.err # error file
#PBS -o ../LogFiles/EWAS_LM.log # output file


## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


cd $PBS_O_WORKDIR

module load R/3.6.0-foss-2019a

Rscript EWAS/Simulations.r rmdConfig.mrc 
