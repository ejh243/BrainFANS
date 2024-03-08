#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=DNAm/logFiles/%u/QTLs-%A_%a.err # error file
#SBATCH --output=DNAm/logFiles/%u/QTLs-%A_%a.log # output file
#SBATCH --job-name=QTLs-%A_%a
#SBATCH --array=1-22


## print start date and time
echo Job started on:
date -u
JOBNAME="QCDNAdata"
    
## needs to be executed from the scripts folder
echo "Changing Folder to: " $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

## format paths in config file with project name
echo "Loading config file for project: " $1
export PROJECT=$1
POP=$2
source ./DNAm/config/config.txt 

module load R/3.6.3-foss-2020a

Rscript DNAm/analysis/QTLs/runQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/ ${SLURM_ARRAY_TASK_ID} covariates_sox10.txt

Rscript DNAm/analysis/QTLs/runCTQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/ ${SLURM_ARRAY_TASK_ID} covariates_sox10.txt

Rscript DNAm/analysis/QTLs/runCTQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/ ${SLURM_ARRAY_TASK_ID} covariates_neun.txt

Rscript DNAm/analysis/QTLs/runQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/NeuN+ ${SLURM_ARRAY_TASK_ID} covariates.txt

Rscript DNAm/analysis/QTLs/runQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/Sox10+ ${SLURM_ARRAY_TASK_ID} covariates.txt

Rscript DNAm/analysis/QTLs/runQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/Double- ${SLURM_ARRAY_TASK_ID} covariates.txt

Rscript DNAm/analysis/QTLs/runQTLByChr.r ${DATADIR}/4_analysis/QTLs/Input/$POP/Total ${SLURM_ARRAY_TASK_ID} covariates.txt