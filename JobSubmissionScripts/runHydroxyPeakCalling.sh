#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=24:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/5hmCPeakCalling.err # error file
#PBS -o LogFiles/5hmCPeakCalling.log # output file

## Output some useful job information

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


source hydroxy/CGEX/config.txt

## index bamfiles
BAMFILES=($(ls ${ALIGNEDDIR}/*_L00.bml.GRCh38.karyo.deduplicated.bam))
module load SAMtools
for f in ${BAMFILES[@]}
do
	samtools index ${f}
done

module load MACS2

cd ${SCRIPTDIR}/hydroxy/CGEX/
./macsPeakCalling.sh

