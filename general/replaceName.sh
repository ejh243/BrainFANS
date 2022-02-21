#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=replaceName-%A.o
#SBATCH --error=replaceName-%A.e
#SBATCH --job-name=replaceName-%A

## this script should be run from the script directory as follows: 
## sbatch path/to/scripts 

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

FILES=()
FILES+=$( find $1 -name *)

echo

python ${SLURM_SUBMIT_DIR}/replaceName.py $2 $3 "${FILES[@]}"

