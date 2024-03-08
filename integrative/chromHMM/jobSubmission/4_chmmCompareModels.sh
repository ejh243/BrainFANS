#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/compareModel-%A.o
#SBATCH --error=integrative/chromHMM/logFiles/compareModel-%A.e
#SBATCH --job-name=compareModel

##this script expects to be submitted from <git-repo> in the format
## sbatch path/to/script <project-name> <refmodel i.e emissions_20.txt>

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
echo 'Loading config file for project: ' $1
INTPROJECT=$1
source ./integrative/chromHMM/config/config.txt

echo "Reference model is: " $2
echo "Output written to: " "${MODELDIR}/compare"

mkdir -p ${MODELDIR}/compare

module load Java

java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar CompareModels ${MODELDIR}/$2 ${MODELDIR} ${MODELDIR}/compare/comparedModel

if [[ $? == 0 ]]
then
	echo "Finished comparison"
fi

#-color 118, 90, 168