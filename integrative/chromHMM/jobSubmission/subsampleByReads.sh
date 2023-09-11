#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=240:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/subsample-%A_%a.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/subsample-%A_%a.e
#SBATCH --job-name=subsample

## print start date and time
echo Job started on:
date -u

start=`date +%s`

## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT

PROJECT=$2

DATA=(${PROJECT//// }) ##split project by / so that data is the datatype

if [[ ${DATA} == 'ChIPSeq' ]]
then 
	DATA=$3
	mark=$(echo $DATA | awk '{print tolower($0)}')
else
	mark="."
fi

source ./integrative/chromHMM/config/config.txt

mkdir -p $INDIR
mkdir -p $OUTDIR
echo $COMPDIR
echo
#-----------------------------------------------------------------------#

cd ${ALIGNEDDIR}
module load SAMtools
module load R/3.6.3-foss-2020a
module load Java/11.0.2
module load BEDTools 

FILES=($(grep "$mark" ${METADIR}/cellMarkFileTable.txt | awk '{print $3}'))
CTLFILES=($(grep "$mark" ${METADIR}/cellMarkFileTable.txt | awk '{print $4}'))

mark=${mark/	/_} #change to e.g. N+_h3k27ac

if [[ ${DATA} == 'ATACSeq' ]] ## use the original bam files (shift later)
then
	FILES=(${FILES[@]%.tn5.tagAlign.gz}) # remove bed endings 
	FILES=(${FILES[@]/%/.filt.nodup.bam}) # add bam ending
fi

infile=downsampled.0.bam
ctrlfile=downsampledCtrl.0.bam

## create the full merged file if it doesn't already exist
if ! test -f ${INDIR}/downsampled.0.bam
then
	echo "Creating joined file"
	samtools merge ${INDIR}/downsampled.0.bam ${FILES[@]}

	if [[ ${#CTLFILES[@]} != 0 ]]; #if ctrlfile length not 0, shouldnt run for atac
	then 
		echo "Creating control file"
		samtools merge ${INDIR}/downsampledCtrl.0.bam ${CTLFILES[@]}
	fi

else
	echo 'Merged file found so skipping merge.'
fi

if ! test -f $OUTDIR/lines_1.txt
then
	totalLines=$(samtools view -c ${INDIR}/${infile} )
	echo -e "0\t"$totalLines >> ${OUTDIR}/lines_${SLURM_ARRAY_TASK_ID}.txt 
else
	totalLines=$(head -n 1 $OUTDIR/lines_1.txt | awk -F'\t' '{print $2}') # catch 
fi	

FILES=(${FILES[@]%.filt.nodup.bam}) #remove chipseq endings
FILES=(${FILES[@]%.nodup.bam}) # remove wgbs endings
FILES=(${FILES[@]%.deduplicated.bam}) # remove hydroxy endings

#avLines=($( Rscript ${SCRIPTDIR}/general/processing/indexVector.R ${PROJECT} ${FILES[@]}))

cd ${INDIR}
#totalLines=$(samtools view -c ${INDIR}/${infile} )

echo "Total number of reads in file:" ${totalLines}

avLines=50000000
echo "Number of reads to downsample by:" ${avLines}


## calculate fraction for the remainder of average reads in total (so that smallest fraction is at the end and not beginning)
fileNo=$(($totalLines / $avLines))
FRACTION=$(bc -l <<< $avLines*$fileNo/$totalLines)

# start downsampling 
for x in $(eval echo "{1..$fileNo}")
do
	echo "Starting downsample on" $infile
	echo 'Fraction of average to whole is:' $FRACTION

	java -jar $EBROOTPICARD/picard.jar DownsampleSam \
	   -I ${infile} \
	   -O ${SLURM_ARRAY_TASK_ID}.downsampled.${x}.bam \
	   -STRATEGY ConstantMemory \
	   -P ${FRACTION} \
       -ACCURACY 0.0001 \
       -VALIDATION_STRINGENCY SILENT

	totalLines=$(samtools view -c ${SLURM_ARRAY_TASK_ID}.downsampled.${x}.bam )
	echo 'New totalLines is:' $totalLines
	echo -e $x'\t'$totalLines >> ${OUTDIR}/lines_${SLURM_ARRAY_TASK_ID}.txt

	infile=${SLURM_ARRAY_TASK_ID}.downsampled.${x}.bam
	echo "File is:" ${infile}
	FRACTION=$(bc -l <<< 1-$avLines/$totalLines)
	echo
done


end=`date +%s`

runtime=$((end-start))

echo "Total runtime is:" ${runtime}

## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}