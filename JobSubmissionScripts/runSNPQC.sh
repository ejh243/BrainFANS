#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=10:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/QCSNPdata.err # error file
#PBS -o LogFiles/QCSNPdata.log # output file


## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


cd $PBS_O_WORKDIR

####### 

## NOTE: Do not store confidenial information in this file use the config file

######



source ./SNPdata/config.txt

#module load PLINK/2.00-alpha1-x86_64
#module load PLINK/1.07-x86_64
module load R/3.5.1-foss-2018b-Python-2.7.15

sh SNPdata/QC.sh
sh SNPdata/CheckEthnicity.sh
# plot PCs
Rscript SNPdata/plotEthnicity.r ${DATADIR}/SNPdata/
sh SNPdata/CheckRelatedness.sh
Rscript SNPdata/plotKinshipCoeff.r ${DATADIR}/SNPdata/

## print finish date and time
echo Job finished on:
date -u
