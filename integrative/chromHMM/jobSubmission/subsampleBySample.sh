#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=240:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/completeness-%A_%a.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/completeness-%A_%a.e
#SBATCH --job-name=completeness

## !! reliant on ${METADIR}/cellMarkFileTable.txt containing the relevant samples

## print start date and time
echo Job started on:
date -u

start=`date +%s`

## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT

PROJECT=$2

if [[ ${PROJECT} == 'ChIPSeq/epiGaba' ]]
then 
	mark="h3k27me3" ## needs to be edited for specificity of the grep search
	DATA="H3K27me3"
else
	mark="."
	DATA=(${PROJECT//// }) ##split project by / so that data is the datatype
fi

source ./integrative/chromHMM/config/config.txt

mkdir -p $COMPDIR
echo $COMPDIR
echo
#-----------------------------------------------------------------------#

cd ${ALIGNEDDIR}
module load R/3.6.3-foss-2020a
module load Java/11.0.2

FILES=($(grep "$mark" ${METADIR}/cellMarkFileTable.txt | awk '{print $3}'))
CTLFILES=($(grep "$mark" ${METADIR}/cellMarkFileTable.txt | awk '{print $4}'))

mark=${mark/	/_} #change to e.g. N+_h3k27ac

# ensure file is empty before appending
echo -n '' > ${METADIR}/cellMarkFileTable_${mark}.txt
echo -n '' > ${COMPDIR}/completeness_${mark}.txt


cd ${ALIGNEDDIR}
fileNo=${#FILES[@]}
echo "Number of files is" $fileNo
fileNo=$((${#FILES[@]} - 1))

for x in $(eval echo "{0..$fileNo}")
do
	## PART 1: BINARISE AND COUNT ZEROS
	echo "Starting binarisation"

	echo -e 'N+\tmark\t'${FILES[x]}'\t'${CTLFILES[x]} >> ${METADIR}/cellMarkFileTable_${mark}.txt

	echo 'Samples are:'
	more ${METADIR}/cellMarkFileTable_${mark}.txt

	echo "Run number is:" $x >> ${COMPDIR}/chromhmm_${mark}.log

	if [[ $DATA == 'ATACSeq' ]]
	then
		java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBed -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable_${mark}.txt ${COMPDIR} &>> ${COMPDIR}/chromhmm_${mark}.log
	else
		java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBam -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable_${mark}.txt ${COMPDIR} &>> ${COMPDIR}/chromhmm_${mark}.log
	fi

	if [[ $? != 0 ]]
	then
		echo "Binarisation incomplete, exiting."
		exit 1
	fi

	echo "Counting marks found"
	# r script which counts number of zeros and writes to a tsv completeness_${mark}.txt
	Rscript ${SCRIPTDIR}/integrative/chromHMM/processing/chmm1count.R ${INTPROJECT} ${DATA} ${mark}

	# remove remaining files
	cd $COMPDIR
	binFILES=$( ls *binary* )
	rm $binFILES
done

end=`date +%s`

runtime=$((end-start))

echo "Total runtime is:" ${runtime}
