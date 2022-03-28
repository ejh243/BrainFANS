#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=100:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=general/logFiles/downloadPsych-%A-%a.o
#SBATCH --error=general/logFiles/downloadPsych-%A-%a.e
#SBATCH --job-name=downloadPsych

#needs to be submitted with the arguments 
#sbatch <path/to/script> <project-including-data-type> i.e. WGBS/epiGaba and <syn-password>

## print start date and time
echo Job started on:
date -u

export PROJECT=$1
echo 'Loading config file for project: ' $PROJECT
source ./sequencing/BSSeq/config/config.txt

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#--------------------------------------------------------------------#

synCode=($(head -n ${SLURM_ARRAY_TASK_ID} ${METADIR}/synCodes.txt | tail -1))

module load R/3.6.0-foss-2019a

Rscript general/processing/downloadPsychEncode.r $2 ${synCode} ${RAWDATADIR}
