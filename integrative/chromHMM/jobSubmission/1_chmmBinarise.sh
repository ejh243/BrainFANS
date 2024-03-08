#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/binarise-%A.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/binarise-%A.e
#SBATCH --job-name=binarise

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

## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT
PROJECT=$2
fileType=$3
tissue=$4

#change project into folder name for binarised folder 
DIR=$(echo $PROJECT | tr / _)

echo "Loading config file for input data project:" $PROJECT
source ./integrative/chromHMM/config/config.txt

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

#check tissue specified and assign default if not
if [[ $4 == '' ]]
then
	echo 'Tissue not specified, using default of prefrontal cortex|PFC'
	tissue="prefrontal cortex|PFC"
fi
echo

#-----------------------------------------------------------------------#

cd $SLURM_SUBMIT_DIR

echo
echo 'Starting creating cellMarkFileTable for' $PROJECT 'at:'
date -u

module purge
module load R/3.6.3-foss-2020a

Rscript ./integrative/chromHMM/processing/chmmMakeCellMarkFile.r ${PROJECT} ${tissue}

if [[ $? == 0 ]]
	then echo "cellMarkFileTable created"

	## start running the binarisation
	echo "Changing to data directory: " ${ALIGNEDDIR}
	cd ${ALIGNEDDIR}

	module purge
	module load Java

	mkdir -p ${BINARISEDIR}/${DIR}

	if [[ $fileType == 'BED' ]];
	then
		echo "Running binarisation on bed files; outdir:" ${BINARISEDIR}/${DIR}

		java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBed -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable.txt ${BINARISEDIR}/${DIR}

		echo "Finished binarisation"
	fi

	if [[ $fileType == 'BAM' ]];
	then
		echo "Running binarisation on bam files; outdir:" ${BINARISEDIR}/${DIR}

		java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBam -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable.txt ${BINARISEDIR}/${DIR}

		echo "Finished binarisation"
	fi
fi


