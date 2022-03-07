#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/logFiles/binarise-%A_%a.o
#SBATCH --error=integrative/logFiles/binarise-%A_%a.e
#SBATCH --job-name=binarise-%A_%a

# This performs binarisation on aligned (and shifted for ATAC) bam files (pre peak calling)
# It expects a specified cellmarkfiletable.txt file directing the call to the files and a method denoting 
# whether selected files are in bed or bam format

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1
fileType=$2

## check file type is specified and matches required options
if [ -z "$fileType" ];
    then 
        { echo "No file type specified using default of bam files" ; fileType="BAM"; }
    else 
        if [ "$fileType" != "BED" ] && [ "$fileType" != "BAM" ];
            then 
                { echo "Unknown file type specified" ; exit 1; }            
    	fi
fi

echo 'Changing to data directory: ' ${ALIGNEDDIR}
cd ${ALIGNEDDIR}

module load Java

if [[ $fileType == 'BED' ]];
then
	echo "Running binarisation on bed files"

	java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBed -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable.txt ${BINARISEDIR}

	echo "Finished binarisation"
fi

if [[ $fileType == 'BAM' ]];
then
	echo "Running binarisation on bam files"

	java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBam -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable.txt ${BINARISEDIR}

	echo "Finished binarisation"
fi
