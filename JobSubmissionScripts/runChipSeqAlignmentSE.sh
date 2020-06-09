#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=150:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/ChipAlignmentAD.err # error file
#PBS -o LogFiles/ChipAlignmentAD.log # output file

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
#cd ${SCRIPTDIR}
source ./ChipSeq/config.ec.ad

## run fastqc
module load FastQC
mkdir -p ${DATADIR}/fastqc
fastqc ${DATADIR}/*q.gz --outdir=${DATADIR}/fastqc

## rn fastp
module load fastp
cd ${DATADIR}
mkdir -p 11_trimmed
mkdir -p 11_trimmed/fastp_reports/
FQFiles=($(ls *q.gz))
for f in ${FQFiles[@]};	
do	
	## find both paired files
	sampleName=${f/_R*}
	if [ ! -f 11_trimmed/fastp_reports/${sampleName}_fastp.json ]	
	then
		outf=${f/.f/_trimmed.f}
		## as SE data no trimming is performed as may interfer with deduplication
		fastp --in1=${f} --out1=11_trimmed/${outf} --length_required=27 --thread=8 --html=11_trimmed/fastp_reports/${sampleName}_fastp.html --json=11_trimmed/fastp_reports/${sampleName}_fastp.json
	fi
done


## merge QC output
module purge ## had conflict issues if this wasn't run first
module load MultiQC/1.2-intel-2017b-Python-2.7.14
cd ${DATADIR}
multiqc . -f ## can add flag to ignore certain folders if needed

## run alignment

module purge ## had conflict issues if this wasn't run first
module load Bowtie2
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
module load BEDTools
#module load Java
cd ${SCRIPTDIR}/ChipSeq/
./alignmentSE.sh ## by using ./ rather than sh executes script in current session and can make use of variables alredy declared.


module purge
module load BEDTools
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
cd ${SCRIPTDIR}/ChipSeq
./calcENCODEQCMetricsSE.sh

echo Starting peak calling at:
date -u
module purge
module load MACS2
cd ${SCRIPTDIR}/ChipSeq/
./peakCallingSE.sh

## print finish date and time
echo Job finished on:
date -u
