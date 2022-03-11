#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=100:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e ../LogFiles/DownloadPsychENCODE.err # error file
#PBS -o ../LogFiles/DownloadPsychENCODE.log # output file

module load R/3.6.0-foss-2019a

cd $PBS_O_WORKDIR

Rscript downloadPsychEncode.r $1
