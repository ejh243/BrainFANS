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

echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME

## print start date and time
echo Job started on:
date -u


cd $PBS_O_WORKDIR

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
