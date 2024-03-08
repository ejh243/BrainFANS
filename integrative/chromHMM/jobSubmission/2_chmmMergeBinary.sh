#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/mergeBinary-%A.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/mergeBinary-%A.e
#SBATCH --job-name=mergeBinary


# This performs mergebinary on the binarised data output
# sbatch /path/to/jobscript/ <project-name>

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder

echo 'Loading config file for project' $1
INTPROJECT=$1
source ./integrative/chromHMM/config/config.txt

mkdir -p ${MERGEDIR}

echo "Running MergeBinary"
module load Java

java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar MergeBinary ${BINARISEDIR} ${MERGEDIR}

if [[ $? == 0 ]]
then
	echo "Finished merge"
fi