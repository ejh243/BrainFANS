#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=general/logFiles/%u/replaceName-%A.o
#SBATCH --error=general/logFiles/%u/replaceName-%A.e
#SBATCH --job-name=replaceName-%A

## this script should be run from the script directory as follows: 
## sbatch replaceName.sh <path/to/scripts> <old-filename> <new-filename>

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder


FILES=()
FILES+=$( find $1 -name '*.*')

echo "Replacing" $2 "with" $3
echo
python ${SLURM_SUBMIT_DIR}/general/processing/replaceName.py $2 $3 "${FILES[@]}"

