#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e QCSNPdata.err # error file
#PBS -o QCSNPdata.log # output file

## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######

cd $PBS_O_WORKDIR

source ./SNPdata/config.txt

module load Perl/5.26.1-foss-2018a

cd ${DATADIR}/scripts ## all scripts are written to be executed from the scripts folder

sh SNPdata/PennCNV.sh

## print finish date and time
echo Job finished on:
date -u
