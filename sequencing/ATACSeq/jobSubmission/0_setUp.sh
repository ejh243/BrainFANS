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
## EXECUTION: sbatch ./sequencing/ATACSeq/jobSubmission/0_setUp.sh <project name>                            ||
## - execute from scripts directory                                                                          ||
##                                                                                                           ||
## INPUTS:                                                                                                   || 
## $1 -> <project name> name of project to set up pipeline for.                                              || 
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


LOG_DIR=ATACSeq/logFiles/${USER}/${SLURM_JOB_ID}
echo "Log files will be moved to dir: " $LOG_DIR
mkdir -p $LOG_DIR
mv ATACSetUpS0-${SLURM_JOB_ID}* $LOG_DIR

## ============ ##
##    CHECKS    ##
## ============ ##

## CONFIG FILE ##

## Check if config file for input project exist

if [[ -s "/lustre/projects/Research_Project-MRC190311/ATACSeq/$1/config.txt" ]]
then
  source "/lustre/projects/Research_Project-MRC190311/ATACSeq/$1/config.txt"
  echo "Config file for $1 project found."
else 
  { echo "Config file for chosen project does not exist or is in the wrong directory. Please create this first." ; exit 1; }  
fi

## Check if R config file for input project exist

if [[ -s "/lustre/projects/Research_Project-MRC190311/ATACSeq/$1/config.r" ]]
then
  echo "R config file for $1 project found."
else 
  { echo "R config file for chosen project does not exist or is in the wrong directory. Please create this first." ; exit 1; }  
fi

## DIRECTORIES ##

## Check that the metadata directory exists and contains sample list.
if [[ ! -d ${METADIR} ]]
then
	echo 'METADATA directory does not exist'
	exit 1
elif  test -f ${METADIR}/samples.txt; then
	mapfile -t SAMPLEIDS < ${METADIR}/samples.txt 
	echo "Sample list found. Number of sample IDs found:" """${#SAMPLEIDS[@]}"""
else
	echo 'samples.txt file does not exist in METADIR'
fi

## Check that the sampleSheet exist
if [[ -d ${METADIR} ]] && [[ -s "${METADIR}/sampleSheet.csv" ]]
then
  SAMPLEIDS2=( $(tail -n +2 ${METADIR}/sampleSheet.csv  | cut -d ',' -f2) )
  echo "Sample sheet found. Number of sample IDs found:"" ""${#SAMPLEIDS2[@]}"""
else
  echo "sampleSheet.csv file does not exist in METADIR."
fi

## Check if remaining main directories exist and create if not
dirs=(${RAWDATADIR} ${FASTQCDIR} ${TRIMDIR} ${ALIGNEDDIR} ${PEAKDIR})
type=("DATA" "FASTQC" "TRIMMED" "ALIGNED" "PEAK CALLED")

for (( x=0; x<${#dirs[@]}; x++ ))
do
	if [ -d ${dirs[x]} ] 
	then 
		echo ${type[x]} "directory:" ${dirs[x]}
	elif [ ${dirs[x]} = $RAWDATADIR ]; then 
		echo ${dirs[x]} "directory does not exist"
		exit 1
	else
		echo "Creating directory: " ${dirs[x]}
		mkdir -p ${dirs[x]} 
	fi
done

## MODULES ##

## Check that the SAMstats environment exists and create it if not

module load Anaconda3
if [ -d ${ENVDIR} ]
then 
	echo "ENCODEQC environment exists"
  source ${CONDAENV}
  conda activate ${ENVDIR}
  echo "List of packages included in conda environment is output in project directory"
  conda list > "/lustre/projects/Research_Project-MRC190311/ATACSeq/$1/condaEnv.txt"
  conda deactivate
else
	echo "ENCODEQC environment does not exist"
	echo "Creating environment"
  source ${CONDAENV}
	conda create --prefix ${ENVDIR}
	conda activate ${ENVDIR}
	echo "Installing SAMstats"
	conda install -c bioconda samstats
	conda deactivate
	echo "ENCODEQC environment has been created and SAMstats installed"
fi

## Check that pip environment exist and create it if not
module load Python/3.9.6-GCCcore-11.2.0-bare

if [[ -d ${PIPENV} ]]
then
  echo "Pip environment exists."
  source ${PIPENV}/bin/activate
  echo "List of packages included in pip environment is output in project directory"
  pip list > "/lustre/projects/Research_Project-MRC190311/ATACSeq/$1/pipEnv.txt"
  deactivate
else
  echo "Pip environment does not exist."
  echo "Creating environment."
  python -m venv ${PIPENV}
  source ${PIPENV}/bin/activate
  echo "Installing basic libraries needed in pipeline"
  pip install Cython
  pip install cykhash
  pip install numpy
  pip install hmmlearn
  pip install scikit-learn
  pip install scipy
  pip install macs3
  deactivate
fi

module purge

echo "R v4.2 is required and will be loaded if found"
## Check that required R version 4.2 exist
module load R/4.2.1-foss-2022a

echo "R packages will be checked and installed if not found."
## Check that required R libraries are installed
Rscript ${SCRIPTDIR}/ATACSeq/preprocessing/setUpRenv.r ${USER}

echo "Set up for project $1 complete."
  
  