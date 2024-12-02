#!/bin/bash
#SBATCH --export=ALL
#SBATCH -p mrcq
#SBATCH --time=24:00:00
#SBATCH -A Research_Project-MRC190311
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END
#SBATCH --error=format1KG.err
#SBATCH --output=format1KG.log

## Output some useful job information

echo "SLURM: working directory is $SLURM_SUBMIT_DIR"
echo "SLURM: job identifier is $SLURM_JOB_ID"
echo "SLURM: job name is $SLURM_JOB_NAME"

## print start date and time
echo Job started on:
date -u


cd "$SLURM_SUBMIT_DIR" || exit 1

####### 
## NOTE: Do not store confidential information in this file use the config file
######

module load BCFtools/1.9-foss-2018b

source ./Misc/config.txt


sh Misc/Download1000GData.sh
sh Misc/Format1000GPlinkgr38.sh

## print finish date and time
echo Job finished on:
date -u
