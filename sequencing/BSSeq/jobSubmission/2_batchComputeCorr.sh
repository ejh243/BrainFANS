#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=BSSeq/logFiles/%u/BSSeqCorrelation-%A_%a.o
#SBATCH --error=BSSeq/logFiles/%u/BSSeqCorrelation-%A_%a.e
#SBATCH --job-name=BSSeqCorrelation

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./BSSeq/config/config.txt 

## check directories
echo "Project directory is: " $DATADIR
echo 'Script directory is: ' ${SCRIPTDIR}

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#-----------------------------------------------------------------------#

module purge
module load R

## get column number of tissue 
column=$(awk '$1 == "fraction"{print NR;exit} ' RS="," ${METADIR}/sampleSheet.csv)

echo 'Column number is ' $column

## get tissue for sample
IFS=$'\n' # need to set this as \n rather than default - space, \t and then \n - so that elements are expanded using \n as delimiter
tissue=($(awk -F"," -v col="$column" '(NR>1) {print $col}' ${METADIR}/sampleSheet.csv | sort -u))
unset IFS
echo 'Tissue is' ${tissue[${SLURM_ARRAY_TASK_ID}]}

# finds all samples with this tissue and save to array
SAMPLEIDS=($(grep "${tissue[${SLURM_ARRAY_TASK_ID}]}" ${METADIR}/sampleSheet.csv | awk -F',' '{print $1}'))

# replace spaces and remove other characters
tissue=${tissue[${SLURM_ARRAY_TASK_ID}]// /_}
tissue=${tissue//[()]/}

echo 'Tissue file name is:' ${tissue}.corr.qc

Rscript ${SCRIPTDIR}/BSSeq/preprocessing/computeCorr.r ${PROJECT} ${tissue} ${SAMPLEIDS[@]}

## move log files into a folder
cd ${SCRIPTDIR}/BSSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv BSSeqCorrelation-${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}
