#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ATACSeq/logFiles/ATACAlignment-%A_%a.o
#SBATCH --error=ATACSeq/logFiles/ATACAlignment-%A_%a.e
#SBATCH --job-name=ATACAlignment-%A_%a.e

## print start date and time
echo Job started on:
date -u

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1 

all=$#

if [[ $2 =~ 'config.dev' ]]
then
	echo "Loading development config file:  "
	echo $2
	source ./$2

	step=$3
	all=1

else
	step=$2
fi

echo 'script directory is: ' ${SCRIPTDIR}

if [ ${all} == 1 ] || [[ ${step} =~ 'FASTQC' ]]
then
	echo "this runs fastqc"
fi
