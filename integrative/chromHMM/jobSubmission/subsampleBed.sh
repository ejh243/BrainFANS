#!/bin/bash
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


FILE=($(find $DATADIR/1_raw -name "*$DATA*"))
CTRLFILE=($(find $DATADIR/1_raw -name "*Input*"))

infile=downsampled.0.bed
ctrlfile=downsampledCtrl.0.bed

## create the full merged file if it doesn't already exist
if ! test -f ${INDIR}/downsampled.0.bed
then
	echo "Creating joined file"
	cp $FILE $INDIR/${infile}.gz	
	zcat $INDIR/${infile}.gz > $INDIR/$infile	

	cp ${CTRLFILE} $INDIR/${ctrlfile}.gz
	zcat $INDIR/${ctrlfile}.gz > $INDIR/$ctrlfile
	rm *.gz
else
	echo 'Merged file found so skipping merge.'
fi

if ! test -f $OUTDIR/lines_1.txt
then
	totalLines=$(wc -l < ${INDIR}/${infile} )
	echo -e "0\t"$totalLines >> ${OUTDIR}/lines_1.txt 
else
	totalLines=$(head -n 1 $OUTDIR/lines_1.txt | awk -F'\t' '{print $2}') # catch 
fi	

cd ${INDIR}

#totalLines=$(samtools view -c ${INDIR}/${infile} )
echo "Total number of reads in file:" ${totalLines}

avLines=10000000
echo "Number of reads to downsample by:" ${avLines}


## calculate fraction for the remainder of average reads in total (so that smallest fraction is at the end and not beginning)
fileNo=$(($totalLines / $avLines))
FRACTION=$(bc -l <<< $avLines*$fileNo) 
echo $FRACTION
if [[ $FRACTION == $totalLines ]]
then
	FRACTION=$(bc -l <<< $totalLines-$avLines)
	fileNo=$(($fileNo-1))
fi

# start downsampling 
for x in $(eval echo "{1..$fileNo}")
do
	echo "Starting downsample on" $infile

	outfile=1.downsampled.${x}.bed

	shuf -n $FRACTION $infile > $outfile

	totalLines=$(wc -l < $outfile)
	echo 'New totalLines is:' $totalLines
	echo -e $x'\t'$totalLines >> ${OUTDIR}/lines_1.txt

	infile=1.downsampled.${x}.bed
	echo "File is:" ${infile}
	FRACTION=$(bc -l <<< $totalLines-$avLines)
	echo 'Number of reads to keep is:' $FRACTION
	echo
done


end=`date +%s`

runtime=$((end-start))

echo "Total runtime is:" ${runtime}

## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}