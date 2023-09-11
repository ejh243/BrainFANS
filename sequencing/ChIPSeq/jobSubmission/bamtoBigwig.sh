#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ChIPSeq/logFiles/%u/bamtoBigWig-%A_%a.o
#SBATCH --error=ChIPSeq/logFiles/%u/bamtoBigWig-%A_%a.e
#SBATCH --job-name=bamtoBigWig

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1

source ./ChIPSeq/config/config.txt 
echo "Project directory is: " $DATADIR

#-----------------------------------------------------------------------#

echo "Changing Folder to Data directory "
echo ${ALIGNEDDIR}

cd ${ALIGNEDDIR}

## get sampleName
sampleName=($(head -n `expr ${SLURM_ARRAY_TASK_ID} + 1` ${METADIR}/samples.txt | tail -1))

## find control to use
column=$(awk -v RS=',' '/controlID/{print NR; exit}' $METADIR/sampleSheet.csv)
control=$(awk -F"," -v col="$column" /$sampleName/'{print $col}' ${METADIR}/sampleSheet.csv)
echo 'Control file specified:' $control

module load deepTools

mkdir -p ${ALIGNEDDIR}/visualisation/

echo
echo "Starting bamCompare on " $sampleName "at: "
date -u
echo 'Control name is: ' $control

bamCompare -b1 ${sampleName}*.nodup.bam \
		-b2 ${control}*.nodup.bam \
		-o ${ALIGNEDDIR}/visualisation/${sampleName}.bw \
		--binSize 20 \
		--normalizeUsing BPM \
		--smoothLength 60 \
		--extendReads 150 \
		--centerReads \
		--scaleFactorsMethod None \
		-p 6 2> ${ALIGNEDDIR}/visualisation/${sampleName}.bamCompare.log


echo "EXITCODE: " $?

## move log files into a folder
cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}
