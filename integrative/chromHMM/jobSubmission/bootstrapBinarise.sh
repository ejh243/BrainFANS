#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/bootBinarise-%A_%a.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/bootBinarise-%A_%a.e
#SBATCH --job-name=bootBinarise

##submit array in increments of two 

## print start date and time
echo Job started on:
date -u

start=`date +%s`


## check chrmm project and input data project
INTPROJECT=$1
echo "ChromHMM project is: " $INTPROJECT



ATAC="ATACSeq/rizzardi"
DNAH="DNAhydroxy/MRC"
CHIP="ChIPSeq/epiGaba"
WGBS="WGBS/rizzardi"

dataTypes=($CHIP $ATAC $DNAH $CHIP $WGBS)

itNo=${SLURM_ARRAY_TASK_ID} 

#-----------------------------------------------------------------------#

source ./integrative/chromHMM/config/config.txt

if [ $# == 1 ] || [[ $2 =~ 'BINARISE' ]]
then
	## for loop to binarise each data type separately. generates the separate datatype dirs in the binarised dir. 
	for PROJECT in ${dataTypes[@]}
	do
		#change project into folder name for binarised folder 
		DIR=$(echo $PROJECT | tr / _)

        cd $SLURM_SUBMIT_DIR
		## rerun source so metadata folder is updated with new project
		source ./integrative/chromHMM/config/config.txt

		echo "Metadata folder is:" ${METADIR}

		sampleNo=$((($(grep 'N+' ${METADIR}/samplesForGroupAnalysis.txt | wc -l ) + 1) / 2))

		echo "Number of samples is:" $sampleNo

		grep 'N+' ${METADIR}/cellMarkFileTable.txt | shuf | split - ${METADIR}/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}_ -l $sampleNo -a 1 --additional-suffix=.txt
		splitIndex=(a b) # create array of the split file identifiers

		cd ${SLURM_SUBMIT_DIR}

		for x in {0..1} # iterate both halves of data split
		do
			itNo=$( expr ${SLURM_ARRAY_TASK_ID} + ${x} ) # iteration number
			mkdir -p ${CHROMDIR}/${itNo}
			module purge
			module load Java

			echo "Files are:" $(awk '{print $3}' ${METADIR}/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}_${splitIndex[${x}]}.txt )
			echo $(awk '{print $3}' ${METADIR}/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}_${splitIndex[${x}]}.txt ) >> ${CHROMDIR}/${itNo}/samples.txt


			source ./integrative/chromHMM/config/config.txt
			mkdir -p ${BINARISEDIR}/${DIR}

			if [[ $PROJECT == $ATAC ]];
			then
				echo "Running binarisation on bed files; outdir:" ${BINARISEDIR}/${DIR}

				java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBed -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}_${splitIndex[${x}]}.txt ${BINARISEDIR}/${DIR}
			else
				echo "Running binarisation on bam files; outdir:" ${BINARISEDIR}/${DIR}

				java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar BinarizeBam -b 200 ${CHROMHMM}/CHROMSIZES/hg38.txt ${ALIGNEDDIR} ${METADIR}/cellMarkFileTable_${SLURM_ARRAY_TASK_ID}_${splitIndex[${x}]}.txt ${BINARISEDIR}/${DIR}
			fi
		done

	done
fi

echo 'EXIT CODE: ' $?

end=`date +%s`

runtime=$((end-start))

echo "Total runtime is:" ${runtime}

## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}
