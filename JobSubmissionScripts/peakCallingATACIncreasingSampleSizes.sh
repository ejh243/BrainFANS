#!/bin/sh
#PBS -V # export all environment variables to the batch job.
#PBS -q mrcq # submit to the serial queue
#PBS -l walltime=150:00:00 # Maximum wall time for the job.
#PBS -A Research_Project-MRC190311 # research project to submit under. 
#PBS -l procs=1 # specify number of processors.
#PBS -m e -M e.j.hannon@exeter.ac.uk # email me at job completion
#PBS -e LogFiles/IncSampleSizePeakCalling.err # error file
#PBS -o LogFiles/IncSampleSizePeakCalling.log # output file

## iteratively increase number of samples and redo peak calling
## then subsmaple number of reads (across samples) and redo peak calling

module load MACS2

ALIGNEDDIR=/gpfs/mrc0/projects/Research_Project-MRC190311/ATACSeq/AlignedData
PEAKDIR=/gpfs/mrc0/projects/Research_Project-MRC190311/ATACSeq/CalledPeaks/AllData/ModelSampleSize


TAGFILES=($(ls ${ALIGNEDDIR}/*/*.tn5.tagAlign.gz))

for NSAMPLES in {2..36}
do

	macs2 callpeak -t ${TAGFILES[@]:0:${NSAMPLES}} --outdir ${PEAKDIR} -n ATACSeqPeakCallingNSamples$NSAMPLES -f BED -g 2.9e9  -B --keep-dup all --shift 100 --extsize 200 --nomodel --broad
	
done

for (( NSAMPLES=2; NSAMPLES<=$NMAX; NSAMPLES++ ))
do
 echo Number of samples: $NSAMPLES of $NMAX
done

FRACTIONS=(N S T)

for type in ${FRACTIONS[@]};
do
	TAGFILES=($(ls ${ALIGNEDDIR}/*/*.tn5.tagAlign.gz | grep -E '\-'${type}\|'_'${type}))
	NMAX=${#TAGFILES[@]}
	for (( NSAMPLES=2; NSAMPLES<=$NMAX; NSAMPLES++ ))
	do
		macs2 callpeak -t ${TAGFILES[@]:0:${NSAMPLES}} --outdir ${PEAKDIR} -n ATACSeqPeakCalling_${type}_NSamples$NSAMPLES -f BED -g 2.9e9  -B --keep-dup all --shift 100 --extsize 200 --nomodel --broad
	done
done

