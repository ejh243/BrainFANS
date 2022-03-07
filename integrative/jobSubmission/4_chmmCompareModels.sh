#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p sq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/logFiles/compareModel-%A.o
#SBATCH --error=integrative/logFiles/compareModel-%A.e
#SBATCH --job-name=compareModel

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
echo 'Loading config file'
source ./integrative/config/config.txt

echo "Reference model is: " $1
echo "Output written to: " "${MODELDIR}/compare/$2"

mkdir -p ${MODELDIR}/compare

module load Java

java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar CompareModels  [ -color 118, 90, 168 ]  ${MODELDIR}/$1 ${MODELDIR} ${MODELDIR}/compare/$2

if [[ $? == 0 ]]
then
	echo "Finished comparison"
fi

#-color 118, 90, 168