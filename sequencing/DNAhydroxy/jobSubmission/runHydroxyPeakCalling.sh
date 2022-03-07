#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=DNAhydroxy/logFiles/CEGX5hmCPeakingCalling.err # error file
#SBATCH --output=DNAhydroxy/logFiles/CEGX5hmCPeakingCalling.log # output file
#SBATCH --job-name=CEGX5hmCPeakingCalling


## print start date and time
echo Job started on:
date -u


source $1

module load R/3.6.3-foss-2020a

cd ${SCRIPTDIR}
Rscript DNAhydroxy/preprocessing/5_createSampleListsForPeakCalling.r ${DATADIR} 

## run peak calling with MACS2
module purge
module load MACS2
module load BEDTools

mkdir -p ${PEAKDIR}/MACS2

sampleSets=($(ls ${METADIR}/sampleLists/PeakCallingInputFiles*))

sh DNAhydroxy/preprocessing/6_macsPeakCallingBySampleType.sh ${sampleSets[${SLURM_ARRAY_TASK_ID}]}

## run peak calling with EPIC2
module purge
module load Miniconda2
source activate epic2
mkdir -p ${PEAKDIR}/EPIC2

sh DNAhydroxy/preprocessing/7_epic2PeakCallingBySampleType.sh ${sampleSets[${SLURM_ARRAY_TASK_ID}]}
source deactivate epic2



