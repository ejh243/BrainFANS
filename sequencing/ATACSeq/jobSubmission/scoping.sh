#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=1 # specify number of processors per node
#SBATCH --output=ATACSeq/logFiles/scoping-%A.o
#SBATCH --error=ATACSeq/logFiles/scoping-%A.e


# This script should check paths to directories, number of files to cross-reference with sample number,
# and necessary package installation for batchRunATACAlignment.sh to run.
# It should be submitted from the scripts directory with the config file to be used with the job submission
# script batchRunATACAlignment.sh <config.txt>

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file:" $1
source ./$1 

## Check downstream directories 
echo "Checking directories, number of files and pipeline stage: "

echo "Project name: " $(basename ${DATADIRPE})

dir=(${SCRIPTDIR} ${RAWDATADIR} ${FASTQCDIR} ${FOLDERTRIM} ${ALIGNEDDIR} ${PEAKDIR})
type=("SCRIPT" "DATA" "FASTQC" "TRIMMED" "ALIGNED" "PEAK CALLED")
pattern=('' '*[rR]1*q.gz' '*_fastqc.zip' '*trimmed*q.gz' '*depDup_q30.bam' '*')

## SCRIPTDIR and RAWDATADIR should be checked for existence and correct file number in RAWDATADIR
for x in {0..1}
do 
	if [ ! -d ${dir[x]} ]
	then 
		echo ${dir[x]} "does not exist"
	else
		echo ${type[x]} "directory:" ${dir[x]}
		if [ ${dir[x]} = ${RAWDATADIR} ] 
		then 
			FILES+=($(find ${RAWDATADIR} -name "${pattern[x]}")) ## this command searches for all fq files within
			echo "Number of read_1.fq.gz files found for alignment:"" ""${#FILES[@]}"""
		fi
	fi
done

## Remaining downstream directories should be created if not in existence and number of files reported to assess pipeline stage
for (( x=2; x<${#dir[@]}; x++ ))
do
	if [ -d ${dir[x]} ] 
	then 
		echo ${type[x]} "directory:" ${dir[x]}
		FILES=()		
		FILES+=($(find ${dir[x]} -name "${pattern[x]}"))
		echo "Number of" ${type[x]} "files: " ${#FILES[@]}
	else 
		echo "Creating directory: " ${dir[x]}
		mkdir -p ${dir[x]} 
	fi
done


## Check that the SAMstats environment exists

echo "Checking that SAMstats environment is set up and SAMstats is installed"

module load Anaconda3
if [ -d ${ENVDIR} ]
then 
	echo "ENCODEQC environment exists"
else
	echo "ENCODEQC environment does not exist"
	echo "Creating environment"
	conda create --prefix ${ENVDIR}
	source activate encodeqc
	echo "Installing SAMstats"
	conda install -c bioconda samstats
	source deactivate
	echo "ENCODEQC environment has been created and SAMstats installed"
fi

module unload Anaconda3

## Check R packages are installed

##module load R/3.6.3-foss-2020a

## move log files into a folder
mkdir -p ATACSeq/logFiles/${USER}
mv ATACSeq/logFiles/scoping* ATACSeq/logFiles/${USER}

