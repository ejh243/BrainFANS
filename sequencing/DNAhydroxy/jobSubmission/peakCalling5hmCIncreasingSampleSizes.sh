#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=LogFiles/IncSampleSizePeakCalling.log
#SBATCH --error=LogFiles/IncSampleSizePeakCalling.log
#SBATCH --job-name=IncSampleSizePeakCalling

## iteratively increase number of samples and redo peak calling

module load MACS2

ALIGNEDDIR=/gpfs/mrc0/projects/Research_Project-MRC190311/CEGX/alignedData
PEAKDIR=/gpfs/mrc0/projects/Research_Project-MRC190311/CEGX/calledPeaks/IncreasingSampleSizes

cd ${ALIGNEDDIR}

SAMPLES=(CEG99-708-09 CEG99-708-21 CEG99-708-44 CEG99-708-58 CEG99-708-54 CEG99-708-59)

# find pulldown and control aligned files
PULLDOWNS=()
CONTROLS=()
for sampleName in ${SAMPLES[@]}
do
	CONTROLS+=($(ls ${sampleName}*IC*.GRCh38.karyo.deduplicated.bam))
	PULLDOWNS+=($(ls ${sampleName}*PC*.GRCh38.karyo.deduplicated.bam))
done

for NSAMPLES in {2..6}
do
	macs2 callpeak -t ${PULLDOWNS[@]:0:${NSAMPLES}} -c ${CONTROLS[@]:0:${NSAMPLES}} -f BAM -g hs -n NeuN5hmCPeakCalling_${NSAMPLES}_Samples -B -q 0.01 --outdir ${PEAKDIR}
done

SAMPLES=(CEG99-708-20 CEG99-708-30 CEG99-708-53 CEG99-708-07 CEG99-708-08 CEG99-708-17)

# find pulldown and control aligned files
PULLDOWNS=()
CONTROLS=()
for sampleName in ${SAMPLES[@]}
do
	CONTROLS+=($(ls ${sampleName}*IC*.GRCh38.karyo.deduplicated.bam))
	PULLDOWNS+=($(ls ${sampleName}*PC*.GRCh38.karyo.deduplicated.bam))
done


for NSAMPLES in {2..6}
do
	macs2 callpeak -t ${PULLDOWNS[@]:0:${NSAMPLES}} -c ${CONTROLS[@]:0:${NSAMPLES}} -f BAM -g hs -n Sox5hmCPeakCalling_${NSAMPLES}_Samples -B -q 0.01 --outdir ${PEAKDIR}
done

CONTROLS=($(ls *IC*.GRCh38.karyo.deduplicated.bam))
PULLDOWNS=($(ls *PC*.GRCh38.karyo.deduplicated.bam))

for NSAMPLES in {2..19}
do
	macs2 callpeak -t ${PULLDOWNS[@]:0:${NSAMPLES}} -c ${CONTROLS[@]:0:${NSAMPLES}} -f BAM -g hs -n All5hmCPeakCalling_${NSAMPLES}_Samples -B -q 0.01 --outdir ${PEAKDIR}
done


