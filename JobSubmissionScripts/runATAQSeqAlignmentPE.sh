#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q sq # submit to the serial queue
#PBS -l walltime=72:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/ATAQAlignment.err # error file
#PBS -o LogFiles/ATAQAlignment.log # output file

## needs to be executed from the scripts folder

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
source ./ATACSeq/config.txt

## run fastqc
module load FastQC 
fastqc ${DATADIRPE}/*q.gz

## rn fastp
module load fastp
cd ${DATADIRPE}
mkdir -p 11_trimmed
mkdir -p 11_trimmed/fastp_reports/
R1Files=($(ls *R1*q.gz))
for f in ${R1Files[@]};	
do	
	## find both paired files
	sampleName=${f/_*}
	if [ ! -f 11_trimmed/fastp_reports/${sampleName}_fastp.json ]	
	then
	  pairedFiles=($(ls ${sampleName}*.gz))
	  f1=${pairedFiles[0]}
	  f2=${pairedFiles[1]}
	  outf1=${f1/.f/_trimmed.f}
	  outf2=${f2/.f/_trimmed.f}
	
	  ## trim adapters only do not trim based on quality
	  fastp --detect_adapter_for_pe --length_required=27 --thread=8 --in1=${f1} --in2=${f2} --out1=11_trimmed/${outf1} --out2=11_trimmed/${outf2} --html=11_trimmed/fastp_reports/${sampleName}_fastp.html --json=11_trimmed/fastp_reports/${sampleName}_fastp.json
	fi
done


## merge QC output
module purge ## had conflict issues if this wasn't run first
module load MultiQC/1.2-intel-2017b-Python-2.7.14
cd ${DATADIRPE}
multiqc . -f ## can add flag to ignore certain folders if needed

## run alignment
module purge ## had conflict issues if this wasn't run first
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
#module load Java
./${SCRIPTDIR}/ATACSeq/alignmentPE.sh

module purge ## had conflict issues if this wasn't run first
module load MultiQC/1.2-intel-2017b-Python-2.7.14
cd ${DATADIR}
multiqc . -f 

module purge
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
module load BEDTools
export PATH="$PATH:/gpfs/mrc0/projects/Research_Project-MRC190311/software/atac_dnase_pipelines/utils/"

./${SCRIPTDIR}/ATACSeq/calcENCODEQCmetricsPE.sh

## shift reads for peak calling
./${SCRIPTDIR}/ATACSeq/shiftAlignedReadsPE.sh

echo Starting peak calling at:
date -u
module purge
module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
./${SCRIPTDIR}/ATACSeq/peakCallingPE.sh

## generate QC metrics
module purge
module load R
Rscript CalcChipQCMetrics.r 

## produce QC report


## print finish date and time
echo Job finished on:
date -u
