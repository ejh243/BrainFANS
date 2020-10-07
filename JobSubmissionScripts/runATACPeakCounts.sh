#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=72:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/PeakCounting.err # error file
#PBS -o LogFiles/PeakCounting.log # output file

## needs to be executed from the scripts folder

## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u

module load HTSeq

export ALIGNEDDIR=/gpfs/mrc0/projects/Research_Project-MRC190311/ATACSeq/AlignedData/[SP]*/
export PEAKDIR=CalledPeaks/AllData/FANSPaper

export BAMFILES=$(ls ${ALIGNEDDIR}/*.filt.nmsrt.nodup.bam) 


sh ATACSeq/countsInPeaks.sh