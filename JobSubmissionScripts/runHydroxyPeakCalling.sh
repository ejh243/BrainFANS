#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/CEGX5hmCPeakingCalling.err # error file
#SBATCH --output=LogFiles/CEGX5hmCPeakingCalling.log # output file
#SBATCH --job-name=CEGX5hmCPeakingCalling


## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

source hydroxy/CGEX/config.txt

cd ${SCRIPTDIR}/hydroxy/CGEX/
module load R/3.6.3-foss-2020a

Rscript createSampleListsForPeakCalling.r config.r 

module purge
module load MACS2


./macsPeakCallingBySampleType.sh

