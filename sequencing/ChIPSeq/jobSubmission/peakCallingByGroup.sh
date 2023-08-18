#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=ChIPSeq/logFiles/%u/ChIPPeakCalling-%A_%a.o
#SBATCH --error=ChIPSeq/logFiles/%u/ChIPPeakCalling-%A_%a.e
#SBATCH --job-name=ChIPPeakCalling

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder
echo "Changing folder to: "
echo $SLURM_SUBMIT_DIR
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
cd $SLURM_SUBMIT_DIR

## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export PROJECT=$1


source ./ChIPSeq/config/config.txt 
echo "Project directory is: " $DATADIR

#check tissue specified and assign default if not
if [[ $tissue == '' ]]
then
	echo 'Tissue not specified, using default of prefrontal cortex|PFC'
	tissue="prefrontal cortex|PFC"
fi
echo

#-----------------------------------------------------------------------#

echo 'Creating samplesForGroupAnalysis.txt based on tissue'

module purge
module load R/3.6.3-foss-2020a

Rscript ../general/processing/makeGroupAnalysisFile.r "ChIPSeq/${PROJECT}" "${tissue}"

echo 'Created samplesForGroupPeaks'

if [[ $2 == 'PASS' ]]
then
	echo 'Group is samples that PASSED QC'
	mapfile -t QC1SAMPLES < ${METADIR}/stage1QCSamples.txt
	export GROUP=stage1Pass
	if [[ ${SLURM_ARRAY_TASK_ID} > 0 ]]
	then
		exit 1
	fi
elif [[ $2 == 'FRACTION' ]]
then
	echo 'Groups are cell fractions'
	## get tissue for sample ( we know the column is 2 beacuse we created the file in an earlier script )
	IFS=$'\n' # need to set this as \n rather than default - space, \t and then \n - so that elements are expanded using \n as delimiter
	fraction=($(awk '{print $2}' ${METADIR}/samplesForGroupAnalysis.txt | sort -u))
	unset IFS
	echo 'Fraction is' ${fraction[${SLURM_ARRAY_TASK_ID}]}

	export GROUP="${fraction[${SLURM_ARRAY_TASK_ID}]}"
else
	echo 'Parameter to differentiate groups by not specified, exiting. Please specify PASS or FRACTION as the third argument.'
	exit 1
fi
shift

if [ $# == 1 ] || [[ $2 =~ 'PEAK' ]]
then
	mkdir -p ${PEAKDIR}/subset

	module purge
	module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14
	module load BEDTools

	cd ${SCRIPTDIR}

	marks=($(awk '{print $3}' $METADIR/samplesForGroupAnalysis.txt | sort -u))
	for mark in ${marks[@]}
	do
		sh ./ChIPSeq/preprocessing/groupPeaks.sh ${mark}
	done
fi

if [[ $? == 0 ]]
	then date -u
	echo "Peaks called"
fi

if [ $# == 1 ] || [[ $2 =~ 'FRIP' ]]
then
	cd ${SCRIPTDIR}

	module purge 
	module load BEDTools

	mkdir -p ${PEAKDIR}/QCOutput/subset/

	sh ./ChIPSeq/preprocessing/calcFripGroup.sh ${QC1SAMPLES[@]}
fi

if [[ $? == 0 ]]
	then date -u
	echo "Fraction of reads in peaks calculated"
fi

## move log files into a folder
cd ${SCRIPTDIR}/ChIPSeq/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.* ${SLURM_ARRAY_JOB_ID}
