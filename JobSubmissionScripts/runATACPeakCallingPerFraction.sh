#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/ATACPeakCalling.err # error file
#SBATCH --output=LogFiles/ATACPeakCalling.log # output file
#SBATCH --job-name=ATACPeakCalling

## needs to be executed from the scripts folder


## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

echo Starting peak calling at:
date -u
module purge
module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14

echo ${ALIGNEDDIR}
echo $PEAKDIR

./ATACSeq/peakCallingByFraction.sh ${ALIGNEDDIR} ${PEAKDIR}


## print finish date and time
echo Job finished on:
date -u
