#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH --time=72:00:00 # Maximum wall time for the job.
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=PIPPASetUpS0-%A.log
#SBATCH --error=PIPPASetUpS0-%A.err
#SBATCH --job-name=PIPPASetUpS0

## ==========================================================================================================##
##                                          PIPPA STEP 0: Set up                                             ##
## ==========================================================================================================##
## EXECUTION: sbatch ./jobSubmission/0_setUp.sh <project directory>                                          ||
##                                                                                                           ||
## INPUTS:                                                                                                   || 
## $1 -> <project directory> path to project's directory                                                     || 
##                                                                                                           ||
## DESCRIPTION: This script sets up the required directories, libraries and modules to run the ATAC-seq      || 
## pipeline. Will check the config file for the project and the initial paths required for the first steps of|| 
## the analysis.                                                                                             ||
##                                                                                                           ||
## REQUIRES:                                                                                                 ||
## - File ${METADIR}/samples.txt that lists sample names.                                                    ||
## - File with requirements to set up conda environment specified in ${CONDA_FILE}                           ||
## - File in ${METADIR} with phenotype information specified in ${SAMPLE_SHEET}                              ||
## - Config.txt, config.r files in <project directory>.                                                      ||
## - The following variables specified in config file: META_DIR, MAIN_DIR, LOG_DIR, RAWDATADIR, ALIGNED_DIR, ||
## SAMPLE_SHEET, TRIM_DIR, PEAK_DIR, CONDA_ENV, CONDA, RSCRIPTS_DIR, RLIB, PROJECT, SCRIPTS_DIR, CONDA_FILE  ||
## - intallLibraries.r file in ${RSCRIPTS_DIR} = ./Rscripts directory                                        || 
## ==========================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## print start date and time
echo Job started on:
date -u

LOG_DIR=${SLURM_JOB_ID}
echo "Log files will be in working directory in folder: " $LOG_DIR
mkdir -p $LOG_DIR
mv PIPPASetUpS0-${SLURM_JOB_ID}* $LOG_DIR

cat <<EOF

Running scripts from: $PWD

EOF

## ============ ##
##    CHECKS    ##
## ============ ##

## CONFIG FILE ##

## Check if config file for input project exist

if [[ -s "$1/config.txt" ]]
then
  source "$1/config.txt"
  echo "Config file for project $PROJECT found."
else 
  { echo "Config file for chosen project does not exist or is in the wrong directory. Please create this first." ; exit 1; }  
fi

## Check if R config file for input project exist

if [[ -s "$1/config.r" ]]
then
  echo "R config file for project $PROJECT found. Should be specified in config.txt file as CONFIGR"
else 
  { echo "R config file for chosen project does not exist or is in the wrong directory. Please create this first." ; exit 1; }  
fi

## DIRECTORIES ##

echo "${META_DIR}"
## Check that the metadata directory exists and contains sample list.
if [[ ! -d "${META_DIR}" ]]
then
	echo 'METADATA directory does not exist'
	exit 1
elif  test -f "${SAMPLES_LIST}"; then
	mapfile -t SAMPLEIDS < "${SAMPLES_LIST}"
	echo "Sample list found. Number of sample IDs found:" """${#SAMPLEIDS[@]}"""
else
	echo 'samples.txt file does not exist in METADIR'
fi

## Check that the sampleSheet exist
if [[ -d "${META_DIR}" ]] && [[ -s "${SAMPLE_SHEET}" ]]
then
  SAMPLEIDS=( $(tail -n +2 "${SAMPLE_SHEET}"  | cut -d ',' -f2) )
  echo "Sample sheet found. Number of sample IDs found:"" ""${#SAMPLEIDS[@]}"""
else
  echo "sampleSheet.csv file does not exist in METADIR."
fi

## Check if remaining main directories exist and create if not
dirs=(${RAWDATADIR} ${FASTQCDIR} ${TRIM_DIR} ${ALIGNED_DIR} ${PEAK_DIR})
type=("DATA" "FASTQC" "TRIMMED" "ALIGNED" "PEAKS CALLED")

for (( x=0; x<${#dirs[@]}; x++ ))
do
	if [ -d ${dirs[x]} ] 
	then 
		echo ${type[x]} "directory:" ${dirs[x]}
	elif [ ${dirs[x]} = "${RAWDATA_DIR}" ]; then 
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

if [[ -d "${CONDA_ENV}" ]]
then 
	echo "${PROJECT} Conda environment exists"
else
	echo "${PROJECT} environment does not exist"
	echo "Creating environment"
  conda activate base
  if [[ ! -f "${CONDA_FILE}" ]]; then
    echo "File with packages to install in Conda environment not found: ${CONDA_FILE}"
    exit 1
  fi

  if conda env create --name "${PROJECT}" --file "${CONDA_FILE}"; then
    echo "${PROJECT} Conda environment created successfully"
  else
    echo "${PROJECT} Conda environment creation failed. Please check the error messages."
    exit 1
  fi

  if [[ ! -d "${CONDA_ENV}" ]]; then
    echo "${PROJECT} Conda environment was not found at ${CONDA_ENV} after creation."
    exit 1
  fi

fi

if [[ -f "${RSCRIPTS_DIR}/installLibraries.r" ]]
then 
  if conda activate "${CONDA_ENV}" ; then

    echo "R packages will be checked and installed if not found."
    echo "R libraries will be installed in ${RLIB}"
  
    ## Check that required R libraries are installed
    Rscript "${RSCRIPTS_DIR}/installLibraries.r" "${RLIB}" "${CONFIGR}"
    conda deactivate
  else
    echo "Could not activate Conda environment: ${CONDA_ENV}"
    exit 1
  fi
fi

echo "Set up for  to be run on ${PROJECT} complete."

conda deactivate