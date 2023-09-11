#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/bootstrap-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/bootstrap-%A_%a.e
#SBATCH --job-name=bootstrap

## print start date and time
echo Job started on:
date -u

start=`date +%s`


## check chrmm project and input data project
INTPROJECT=$1
echo "Segway project is: " $INTPROJECT

ATAC="ATACSeq/rizzardi"
DNAH="DNAhydroxy/MRC"
CHIP="ChIPSeq/epiGaba"
WGBS="WGBS/rizzardi"

dataTypes=($ATAC $DNAH $CHIP $WGBS)

export GROUP=N+

#-----------------------------------------------------------------------#

for PROJECT in ${dataTypes[@]}
do
	echo 'Project is' $PROJECT
	if [[ $PROJECT == $CHIP ]]
	then
		source ./integrative/segway/config/config.txt 

		marks=$(awk '{print $3}' $METADIR/samplesForGroupAnalysis.txt | sort -u)

		for mark in ${marks}
		do
			grep ${mark} $METADIR/samplesForGroupAnalysis.txt > $METADIR/samplesForGroupAnalysis_${mark}.txt

			groupFile=${METADIR}/samplesForGroupAnalysis_${mark} #e.g. samplesForGroupAnalysis_H3K27ac.txt
			sampleNo=$((($(grep 'N+' ${groupFile}.txt | wc -l ) + 1) / 2))

			##create two samplesForGroupAnalysis files with half the data each, called e.g. samplesForGroupAnalysis_H3K27ac_1_a.txt
			grep 'N+' ${groupFile}.txt | shuf | split - ${groupFile}_${SLURM_ARRAY_TASK_ID}_ -l ${sampleNo} -a 1 --additional-suffix=.txt
			splitIndex=(a b) # create array of the split file identifiers

			echo 'Groupfile is' ${groupFile}.txt

			for x in {0..1} # iterate both halves of data split
			do
				FILE=${METADIR}/samplesForGroupAnalysis_${mark}_${SLURM_ARRAY_TASK_ID}_${splitIndex[x]}.txt
				echo 'SampleFile is' ${FILE}
				

				itNo=$( expr ${SLURM_ARRAY_TASK_ID} + ${x} ) # iteration number
				mkdir -p ${SEGDIR}/${itNo}

				source ./integrative/segway/config/config.txt

				sh ./sequencing/ChIPSeq/preprocessing/groupPeaks.sh ${FILE} ${mark}

				awk '{print $1,$2,$3,$7}' ${PEAKGROUP}/${GROUP}_${mark}.*.filt > ${PEAKGROUP}/${GROUP}_${mark}.filt.srt.bg
			done
		done
	else
		## rerun source so PROJECT paths are correct
		source ./integrative/segway/config/config.txt 

		sampleNo=$((($(grep 'N+' ${METADIR}/samplesForGroupAnalysis.txt | wc -l ) + 1) / 2))

		grep 'N+' ${METADIR}/samplesForGroupAnalysis.txt | shuf | split - ${METADIR}/samplesForGroupAnalysis_${SLURM_ARRAY_TASK_ID}_ -l ${sampleNo} -a 1 --additional-suffix=.txt
		splitIndex=(a b) # create array of the split file identifiers

		for x in {0..1} # iterate both halves of data split
		do
			itNo=$( expr ${SLURM_ARRAY_TASK_ID} + ${x} ) # iteration number
			mkdir -p ${SEGDIR}/${itNo}

			#rerunning so that itNo is correct
			cd ${SLURM_SUBMIT_DIR}
			source ./integrative/segway/config/config.txt 

			FILE=${METADIR}/samplesForGroupAnalysis_${SLURM_ARRAY_TASK_ID}_${splitIndex[x]}.txt

			module purge
			module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
			module load BEDTools

			if [[ $PROJECT == $ATAC ]]
			then
				sh ./sequencing/ATACSeq/preprocessing/groupPeaks.sh ${FILE}

				## create peakfile with only columns 1,2,3,7 in the format:
				# chromosome [s] chrstart [s] chrstop [s] signalValue (measure of total overall enrichment for the region)
				awk '{print $1,$2,$3,$7}' ${PEAKGROUP}/${GROUP}_atac.*.filt > ${PEAKGROUP}/${GROUP}_atac.filt.srt.bg

			elif [[ $PROJECT == $DNAH ]]
			then
				sh ./sequencing/DNAhydroxy/preprocessing/groupPeaks.sh ${FILE}

				## create peakfile with only columns 1,2,3,7 in the format:
				# chromosome [s] chrstart [s] chrstop [s] signalValue (measure of total overall enrichment for the region)
				awk '{print $1,$2,$3,$7}' ${PEAKGROUP}/${GROUP}_5hmC.*.filt > ${PEAKGROUP}/${GROUP}_5hmC.filt.srt.bg

			elif [[ $PROJECT == $WGBS ]]
			then 
				QC1SAMPLES=($(awk '{print $1}' ${FILE} ))
				FILES=( "${QC1SAMPLES[@]/%/*bismark.cov.gz}" )

				cd ${METHYLDIR} 
				zcat ${FILES[@]} | sort -k1,1 -k2,2n | bedtools map -a ${REF}/genome.window.fa.bg -b - -c 4 -o mean | \
					grep -P "\d$" > ${PEAKGROUP}/'N+_5mC.window.bg'
			fi
		done
	fi
done

## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromHMM/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}