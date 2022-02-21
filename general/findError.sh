#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=findError-%A.o
#SBATCH --error=findError-%A.e
#SBATCH --job-name=findError-%A

## this script should be run from the script directory as follows: 
## sbatch path/to/script <config.txt> <error-phrase-to-search> <path-to-logfiles-from-scriptdir>


## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1
error="$2"
LOGDIR=$( dirname $3 )
logfile=$( basename $3 )*.e

echo "Error phrase is: " "${error}"
echo "Logfile directory is: " $LOGDIR
cd $LOGDIR

#create array of logfiles to pass to script
EFILES=()
EFILES+=($( ls ${logfile} ))
EFILES=( "${EFILES[@]%.e}" ) #remove all file extensions from the files in the list 
echo "Number of error log files found: " ${#EFILES[@]}
echo 'File names are: ' "${EFILES[@]}"

python ${SLURM_SUBMIT_DIR}/findError.py "${error}" "${EFILES[@]}"