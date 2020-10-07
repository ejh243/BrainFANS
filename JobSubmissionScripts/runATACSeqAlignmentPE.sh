#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/ATAQAlignment.err # error file
#SBATCH --output=LogFiles/ATAQAlignment.log # output file
#SBATCH --job-name=ATAQAlignment

## needs to be executed from the scripts folder


## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

echo "Loading config file: "
echo $1
source ./$1 ## config file provided on command line when submitting job

echo "Changing Folder to Data directory "
echo ${DATADIRPE}

cd ${DATADIRPE}
RAWDATADIR=($(find . -name '01_raw_reads'))

## run fastqc
echo "Running FASTQC"
module load FastQC 
for FOLDER in ${RAWDATADIR[@]}
do
	echo "Processing folder: "
	echo ${FOLDER}
	FQFILES=($(find ${FOLDER} -name '*q.gz')) ## this command searches for all fq files within data directory
	fastqc ${FQFILES[@]}
done
date -u
echo "FASTQC Done"


## run fastp if not already run
echo "Running FASTP"
module load fastp
for FOLDER in ${RAWDATADIR[@]}
do
	echo "Processing folder: "

	FOLDERTRIM=${${FOLDER}/#01_raw_reads/11_trimmed}
	echo ${FOLDERTRIM}
	mkdir -p ${FOLDERTRIM}
	mkdir -p ${FOLDERTRIM}/fastp_reports/
	R1Files=($(find ${FOLDER} -name '*[rR]1*q.gz'))
	for f in ${R1Files[@]};	
	do	
		## find both paired files
		sampleName=${f/_*}
		if [ ! -f ${FOLDERTRIM}/fastp_reports/${sampleName}_fastp.json ]	
		then
		  pairedFiles=($(ls ${sampleName}*.gz))
		  f1=${pairedFiles[0]}
		  f2=${pairedFiles[1]}
		  outf1=${f1/.f/_trimmed.f}
		  outf2=${f2/.f/_trimmed.f}
		
		  ## trim adapters only do not trim based on quality
		  fastp --detect_adapter_for_pe --length_required=27 --thread=8 --in1=${f1} --in2=${f2} --out1=${FOLDERTRIM}/${outf1} --out2=${FOLDERTRIM}/${outf2} --html=${FOLDERTRIM}/fastp_reports/${sampleName}_fastp.html --json=${FOLDERTRIM}/fastp_reports/${sampleName}_fastp.json
		fi
	done
done
date -u
echo "FASTP done"

## merge QC output
module purge ## had conflict issues if this wasn't run first
module load MultiQC/1.2-intel-2017b-Python-2.7.14
multiqc . -f ## can add flag to ignore certain folders if needed

## run 
echo "Running alignment"
module purge ## had conflict issues if this wasn't run first
module load Bowtie2/2.3.4.2-foss-2018b
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
#module load Java
cd ${SCRIPTDIR}
./ATACSeq/alignmentPE.sh

date -u
echo "Alignment done"


module purge ## had conflict issues if this wasn't run first
module load MultiQC/1.2-intel-2017b-Python-2.7.14
cd ${ALIGNEDDIR}
multiqc . -f 

module purge
module load SAMtools
module load picard/2.6.0-Java-1.8.0_131
module load BEDTools
export PATH="$PATH:/gpfs/mrc0/projects/Research_Project-MRC190311/software/atac_dnase_pipelines/utils/"

echo "Calculating ENCODE QC metrics"
cd ${SCRIPTDIR}
./ATACSeq/calcENCODEQCMetricsPE.sh
date -u
echo "ENCODE QC metrics calculated"


## shift reads for peak calling
echo "Shifting reads"

cd ${SCRIPTDIR}
./ATACSeq/shiftAlignedReadsPE.sh

date -u
echo "Reads shifted"


echo Starting peak calling at:
date -u
module purge
module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
cd ${SCRIPTDIR}/
./ATACSeq/peakCallingPE.sh


## print finish date and time
echo Job finished on:
date -u
