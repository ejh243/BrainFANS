#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/logFiles/learnModel-%A_%a.o
#SBATCH --error=integrative/logFiles/learnModel-%A_%a.e
#SBATCH --job-name=learnModel

# This performs learnmodel on merged binarisation files 
# to be submitted from <repo> as sbatch --array=1-_ integrative/jobSubmission/3_chmmLearnModel.sh 

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

echo 'Loading config file'
source ./integrative/config/config.txt

if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then
	echo 'Job does not appear to be an array. Please set array number'
	exit 1
else
	echo 'Number of states: ' ${SLURM_ARRAY_TASK_ID}
fi

mkdir -p ${MODELDIR}

module load Java

echo 
echo 'Learning model'
java -mx4000M -jar ${CHROMHMM}/ChromHMM.jar LearnModel -p 0 ${MERGEDIR}/2_output ${MODELDIR} ${SLURM_ARRAY_TASK_ID} hg38

if [[ $? == 0 ]]
then 
	echo "Model(s) learned"
fi

echo 'EXIT CODE: ' $?