#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
# Temporary log files, job ID given to file names in case job is ran multiple
# times in parallel
#SBATCH --error=LDSCNeuralPeaks_%j.err
#SBATCH --output=LDSCNeuralPeaks_%j.log
#SBATCH --job-name=LDSCNeuralPeaks

## ========= ##
##   SETUP   ##
## ========= ##

echo Job started on:
date -u

# This cd is in place so that the script can be ran from anywhere.
# The called scripts in this script are in the same directory as this one.
script_directory=$(dirname "$0")
cd "${script_directory}/.." || exit 1

annotation_flag=false

while getopts "r" OPT; do
  case "$OPT" in
    r )   annotation_flag=true ;;
    \? )  { echo "invalid argument parsed: $OPTARG"; exit 1; } ;;
  esac
done
shift $((OPTIND-1))

configuration_directory=$1


source "${configuration_directory}/config.txt"

mv "${SLURM_SUBMIT_DIR}/LDSCNeuralPeaks_${SLURM_JOB_ID}.log" \
"${LOG_DIR}/LDSCNeuralPeaks_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/LDSCNeuralPeaks_${SLURM_JOB_ID}.err" \
"${LOG_DIR}/LDSCNeuralPeaks_${SLURM_JOB_ID}.err"

## ======== ##
##   MAIN   ##
## ======== ##

if [ $annotation_flag ]; then
  module purge
  module load BEDOPS
  module load BEDTools
  module load R/3.6.3-foss-2020a
  Rscript processing/createAnnotationFiles.r \
  "${configuration_directory}/config.r"
fi

bash jobSubmission/runPartionedHeritabilityOnPeaks.sh

echo Job finished at:
date -u