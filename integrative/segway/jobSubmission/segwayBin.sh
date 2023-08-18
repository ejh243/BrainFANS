#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/window-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/window-%A_%a.e
#SBATCH --job-name=window

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export INTPROJECT=$1
export PROJECT=$2

source ./integrative/segway/config/config.txt 
echo "Project directory is: " $DATADIR

#check tissue specified and assign default if not
if [[ $tissue == '' ]]
then
	echo 'Tissue not specified, using default of prefrontal cortex|PFC'
	tissue="prefrontal cortex|PFC"
fi
echo

#-----------------------------------------------------------------------#

echo 'Creating samplesForGroupAnalysis.txt based on tissue'

#module purge
#module load R/3.6.3-foss-2020a

#Rscript ../general/processing/makeGroupAnalysisFile.r "${PROJECT}" "${tissue}"

fraction=($(awk '{print $2}' ${METADIR}/samplesForGroupAnalysis.txt | sort -u))
echo 'Fraction is' ${fraction[${SLURM_ARRAY_TASK_ID}]}

QC1SAMPLES=($(grep "${fraction[${SLURM_ARRAY_TASK_ID}]}" ${METADIR}/samplesForGroupAnalysis.txt | awk '{print $1}'))

FILES=( "${QC1SAMPLES[@]/%/*bismark.cov.gz}" )


echo 'Files to merge are: ' ${FILES[@]}

module purge
module load BEDTools

echo
echo 'Starting to map file at:'
date -u
 

cd ${METHYLDIR} 

zcat ${FILES[@]} | sort -k1,1 -k2,2n > ${fraction[${SLURM_ARRAY_TASK_ID}]}.bg

bedtools map -a ${REF}/genome.window.fa.bg -b ${fraction[${SLURM_ARRAY_TASK_ID}]}.bg -c 4 -o mean | grep -P "\d$" > ${LOADDIR}/${fraction[${SLURM_ARRAY_TASK_ID}]}'_5mC.window.bg'