#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=1 # specify number of processors per node
#SBATCH --output=ATACSeq/logFiles/%u/scoping-%A.o
#SBATCH --error=ATACSeq/logFiles/%u/scoping-%A.e


# This script should check paths to directories, number of files to cross-reference with sample number,
# and necessary package installation for batchRunATACAlignment.sh to run.
# It should be submitted from the scripts directory with the config file to be used with the job submission
# script batchRunATACAlignment.sh 
# sbatch ../general/jobSubmission/scoping.sh <config.txt>

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder


## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ATACSeq/config/config.txt 
echo "Project directory is: " $DATADIR

## Check downstream directories 
echo "Checking directories, number of files and pipeline stage: "

echo "Project name: " $PROJECT

dir=(${RAWDATADIR} ${FASTQCDIR} ${TRIMDIR} ${ALIGNEDDIR} ${PEAKDIR})
type=("DATA" "FASTQC" "TRIMMED" "ALIGNED" "PEAK CALLED")
pattern=('*q.gz' '*_fastqc.zip' '*q.gz' '*depDup_q30.bam' '*')


## Check that the metadata directory exists and contains the sample list with correct number
if [ ! -d ${METADIR} ]
then
	echo 'METADATA directory does not exist'
	exit 1
elif  test -f ${METADIR}/samples.txt; then
	mapfile -t SAMPLEIDS < ${METADIR}/samples.txt 
	echo "Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
else
	echo 'samples.txt file does not exist in METADIR'
fi


## Remaining downstream directories should be created if not in existence and number of files reported to assess pipeline stage
for (( x=0; x<${#dir[@]}; x++ ))
do
	if [ -d ${dir[x]} ] #if directory exists
	then 
		echo ${type[x]} "directory:" ${dir[x]}
		FILES=()
		for i in ${SAMPLEIDS[@]}
		do 
			FILES+=($(find ${dir[x]} -maxdepth 1 -name "${i}${pattern[x]}"))
			#echo ${FILES[@]}
		done
		echo "Number of" ${type[x]} "files: " ${#FILES[@]}
	elif [ ${dir[x]} = $RAWDATADIR ]; then ##
		echo ${dir[x]} "directory does not exist"
		exit 1
	else 
		echo "Creating directory: " ${dir[x]}
		mkdir -p ${dir[x]} 
	fi
done


## Check that the SAMstats environment exists and create it if not

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

## Check R packages are installed

##module load R/3.6.3-foss-2020a

## move log files into a folder
cd ${SCRIPTDIR}/ATACSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv scoping-${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}/