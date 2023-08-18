#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=400:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=DNAm/logFiles/%u/CTDeconvolutionSimulations-%A_%a.err # error file
#SBATCH --output=DNAm/logFiles/%u/CTDeconvolutionSimulations-%A_%a.log # output file
#SBATCH --job-name=CTModel_%a


# print start date and time
echo Job started on:
date -u
    
# needs to be executed from the scripts folder
echo "Changing Folder to: " $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

# load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1
source ./DNAm/config/config.txt 

module load R/3.6.3-foss-2020a

Rscript DNAm/analysis/neuralCellComposition/fitModelsSimulations.r ${DATADIR}/3_normalised/normalisedMergeWEpigaba.rdata ${REFDIR}/450K_reference DNAm/analysis/neuralCellComposition/NeuralDeconvolutionPanels.csv ${SLURM_ARRAY_TASK_ID}


# print end date and time
echo Job finished on:
date -u
    
