#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSetUpS0-%A.log
#SBATCH --error=ATACSetUpS0-%A.err
#SBATCH --job-name=ATACSetUpS0

## ==========================================================================================================##
##                               ATAC-seq pipeline STEP 0: Set up                                            ##
## ==========================================================================================================##
## EXECUTION: sbatch ./sequencing/ATACSeq/jobSubmission/0_setUp.sh <project directory>                       ||
## - execute from scripts directory                                                                          ||
##                                                                                                           ||
## INPUTS:                                                                                                   || 
## $1 -> <project directory> path to project's directory                                                     || 
##                                                                                                           ||
## DESCRIPTION: This script sets up the required directories, libraries and modules to run the ATAC-seq      || 
## pipeline. Will check the config file for the project and the initial paths required for the first steps of|| 
## the analysis.                                                                                             ||
##                                                                                                           ||
## REQUIRES:                                                                                                 ||
## - File in ${METADIR}/samples.txt that lists sample names.                                                 ||
## ==========================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

LOG_DIR=./ATACSeq/logFiles/${USER}/${SLURM_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACSetUpS0-${SLURM_JOB_ID}* $LOG_DIR

cat <<EOF

Running scripts from: $PWD
Log files will be in current directory: $LOG_DIR

EOF

## ============ ##
##    CHECKS    ##
## ============ ##

## CONFIG FILE ##

## Check if config file for input project exist

if [[ -s "$1/config.txt" ]]
then
  source "$1/config.txt"
  echo "Config file for project $PROJECT  found."
else 
  { echo "Config file for chosen project does not exist or is in the wrong directory. Please create this first." ; exit 1; }  
fi

## Check if R config file for input project exist

if [[ -s "$1/config.r" ]]
then
  echo "R config file for project $PROJECT found."
else 
  { echo "R config file for chosen project does not exist or is in the wrong directory. Please create this first." ; exit 1; }  
fi

## DIRECTORIES ##

echo ${META_DIR}
## Check that the metadata directory exists and contains sample list.
if [[ ! -d ${META_DIR} ]]
then
	echo 'METADATA directory does not exist'
	exit 1
elif  test -f ${META_DIR}/samples.txt; then
	mapfile -t SAMPLEIDS < ${META_DIR}/samples.txt 
	echo "Sample list found. Number of sample IDs found:" """${#SAMPLEIDS[@]}"""
else
	echo 'samples.txt file does not exist in METADIR'
fi

## Check that the sampleSheet exist
if [[ -d ${META_DIR} ]] && [[ -s "${META_DIR}/sampleSheet.csv" ]]
then
  SAMPLEIDS=( $(tail -n +2 ${META_DIR}/sampleSheet.csv  | cut -d ',' -f2) )
  echo "Sample sheet found. Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
else
  echo "sampleSheet.csv file does not exist in METADIR."
fi

## Check if remaining main directories exist and create if not
dirs=(${RAWDATADIR} ${FASTQCDIR} ${TRIM_DIR} ${ALIGNED_DIR} ${PEAK_DIR})
type=("DATA" "FASTQC" "TRIMMED" "ALIGNED" "PEAK CALLED")

for (( x=0; x<${#dirs[@]}; x++ ))
do
	if [ -d ${dirs[x]} ] 
	then 
		echo ${type[x]} "directory:" ${dirs[x]}
	elif [ ${dirs[x]} = ${RAWDATA_DIR} ]; then 
		echo ${dirs[x]} "directory does not exist"
		exit 1
	else
		echo "Creating directory: " ${dirs[x]}
		mkdir -p ${dirs[x]} 
	fi
done

## Activate conda environment with packages/modules
source ${CONDA} 

## Check that a conda environment exists and create it if not

if [[ -d ${CONDA_ENV} ]]
then 
	echo "${PROJECT} conda environment exists"
else
	echo "${PROJECT} environment does not exist"
	echo "Creating environment"
  conda activate base
  if [[ -f ${CONDA_FILE} ]]
  then
	  conda env create --name ${PROJECT} --file ${CONDA_FILE}
    conda deactivate
  else
    echo "File with packages to install in conda environment not found."
  fi
  if [[ -f ${PIP_FILE} ]]
  then
    conda activate ${CONDA_ENV}
    pip install -r ${PIP_FILE}
    conda deactivate
  else
    echo "File with packages to install in conda environment through pip not found."
  fi
  
  if [[ -f ${RSCRIPTS_DIR}/intallLibraries.r ]]
  then 
    conda activate ${CONDA_ENV}
  
cat <<EOF

R packages will be checked and installed if not found.
R libraries will be installed in ${RLIB}

EOF  
    ## Check that required R libraries are installed
    Rscript ${RSCRIPTS_DIR}/intallLibraries.r ${RLIB}
  else
    echo "Setup R libraries script not found. Please make sure this is in the ${RSCRIPTS_DIR} directory"
  fi

  if [[ -d ${CONDA_ENV} ]]
  then 
	  echo "${PROJECT} conda environment has been created"
  else
    echo "${PROJECT} conda environment has not been successfully created. Please check log error file and try again."
  fi
fi

echo "Set up for ATAC-seq pipeline to be run on ${PROJECT} complete."

conda deactivate

