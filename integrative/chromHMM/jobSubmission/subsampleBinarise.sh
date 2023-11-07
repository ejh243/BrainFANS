#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=240:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/subBin-%A_%a.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/subBin-%A_%a.e
#SBATCH --job-name=subBin

## print start date and time
echo Job started on:
date -u

start=`date +%s`

## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT

PROJECT=$2

array=$3. #number at the start of each binarisation

# set array number to nothing if dealing with the first all merged file
if [[ ${SLURM_ARRAY_TASK_ID} == '0' ]]; then array=''; fi

DATA=(${PROJECT//// }) ##split project by / so that data is the datatype

if [[ ${DATA} == 'ChIPSeq' ]]
then 
	DATA=$4
fi

source ./integrative/chromHMM/config/config.txt

mkdir -p $COMPDIR/0_metadata
mkdir -p ${TMPDIR}/${SLURM_ARRAY_TASK_ID}
echo $COMPDIR
echo
#-----------------------------------------------------------------------#

#INDIR=$COMPDIR
#TMPDIR=$COMPDIR

cd ${INDIR}
module purge
module load Java/11.0.2
module load R/3.6.3-foss-2020a


## PART 1: BINARISE AND COUNT ZEROS
echo "Starting binarisation"

echo -n '' > $COMPDIR/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt

##get the filetype to specify cellmark file and also the chromhmm commmand
filetype=$(find -name "downsampled.0.*")
filetype=$(echo ${filetype##*.})

echo 'Filetype is' $filetype

if [[ $filetype == 'bam' ]]; then command='BinarizeBam'; else if [[ $filetype == 'bed' ]]; then command='BinarizeBed'; fi; fi

cd $INDIR
ctrlfile=$(find downsampledCtrl.0.*)
ctrlfile=$(echo '\t'${ctrlfile})

if [[ ${DATA} != 'ATACSeq' ]]
then 
	echo -e 'N+\tmark\t'${array}'downsampled.'${SLURM_ARRAY_TASK_ID}'.'${filetype}${ctrlfile} > $COMPDIR/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt

	echo "Sample is:"
	more $COMPDIR/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt

	java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar $command -b 200 ${CHROMHMM}/CHROMSIZES/hg19.txt ${INDIR} \
						${COMPDIR}/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt ${TMPDIR}/${SLURM_ARRAY_TASK_ID} &>> ${TMPDIR}/${SLURM_ARRAY_TASK_ID}/chromhmm.log
	if [[ $? != 0 ]]
	then
		echo "Binarisation incomplete, exiting."
		exit 1
	fi
else
	echo 'Starting shift'
	SAMPLE=$array'downsampled.'${SLURM_ARRAY_TASK_ID}'.bam'
	sh $SCRIPTDIR/integrative/chromHMM/processing/shiftATACSubsample.sh ${SAMPLE}

	echo -e 'N+\tmark\t'$array'downsampled.'${SLURM_ARRAY_TASK_ID}'.tn5.tagAlign.gz' > $COMPDIR/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt

	echo "Sample is:"
	more $COMPDIR/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt

	java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBed -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${INDIR} \
						${COMPDIR}/0_metadata/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}.txt ${TMPDIR}/${SLURM_ARRAY_TASK_ID} &>> ${TMPDIR}/${SLURM_ARRAY_TASK_ID}/chromhmm.log
	if [[ $? != 0 ]]
	then
		echo "Binarisation incomplete, exiting."

		## move log files into a folder
		cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
		mkdir -p ${SLURM_ARRAY_JOB_ID}
		mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}
		exit 1
	fi
fi

echo "Counting marks found"
# r script which counts number of zeros and writes to a tsv completeness_${SLURM_ARRAY_TASK_ID}.txt
Rscript ${SCRIPTDIR}/integrative/chromHMM/processing/chmm1count.R ${INTPROJECT} ${DATA} ${SLURM_ARRAY_TASK_ID}

# remove remaining files
#cd ${COMPDIR}/${SLURM_ARRAY_TASK_ID}
#binFILES=$( ls *binary* )
#rm $binFILES


cd $COMPDIR

compFiles=($(find . -name completeness.txt | sort -t / -k 2 -g))

if [[ ${#compFiles[@]} == ${SLURM_ARRAY_TASK_COUNT} ]]
then
	echo 'Combining completeness'
	cat ${compFiles[@]} > ${OUTDIR}/completeness_${array}txt
fi

end=`date +%s`

runtime=$((end-start))

echo "Total runtime is:" ${runtime}

## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}