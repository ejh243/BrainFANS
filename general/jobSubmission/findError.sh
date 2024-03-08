#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p pq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=general/logFiles/findError-%A.o
#SBATCH --error=general/logFiles/findError-%A.e
#SBATCH --job-name=findError-%A

## this script should be run from the script directory as follows: 
## sbatch path/to/script <error-phrase-to-search> <path-to-logfiles-from-scriptdir>


## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

error="$1"
LOGDIR=$( dirname $2 )
logfile=$( basename $2 )*.e

echo "Error phrase is: " "${error}"
echo "Logfile directory is: " $LOGDIR
cd $LOGDIR

#create array of logfiles to pass to script
EFILES=()
EFILES+=($( ls ${logfile} ))
EFILES=( "${EFILES[@]%.e}" ) #remove all file extensions from the files in the list 
echo "Number of error log files found: " ${#EFILES[@]}
echo 'File names are: ' "${EFILES[@]}"

python ${SLURM_SUBMIT_DIR}/general/processing/findError.py "${error}" "${EFILES[@]}"