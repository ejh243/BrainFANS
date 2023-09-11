#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/segway/logFiles/%u/gnmdatGenerate-%A_%a.o
#SBATCH --error=integrative/segway/logFiles/%u/gnmdatGenerate-%A_%a.e
#SBATCH --job-name=gnmdatGenerate

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## check chrmm project and input data project
INTPROJECT=$1
echo "Segway project is: " $INTPROJECT
PROJECT=$2


echo "Loading config file"
source ./integrative/segway/config/config.txt
echo

#-----------------------------------------------------------------------#

module load Anaconda3
source activate segway

mkdir -p ${LOADDIR}

DATA=(${PROJECT//// }) ##split project by / so that data is the datatype

if [[ 'WGBS,oxBS' =~ $DATA ]]
then 
	echo "Data has been bisulfite sequenced"
	fraction=($(awk '{print $2}' ${METADIR}/samplesForGroupAnalysis.txt | sort -u))
	echo 'Fraction is' ${fraction[${SLURM_ARRAY_TASK_ID}]}

	QC1SAMPLES=($(grep "${fraction[${SLURM_ARRAY_TASK_ID}]}" ${METADIR}/samplesForGroupAnalysis.txt | awk '{print $1}'))
	FILES=( "${QC1SAMPLES[@]/%/*bismark.cov.gz}" )

	echo 'Files to merge are: ' ${FILES[@]}

	module purge
	module load BEDTools

	echo 'Changing to methylation directory' $METHYLDIR
	cd ${METHYLDIR} 

	zcat ${FILES[@]} | sort -k1,1 -k2,2n | bedtools map -a ${REF}/genome.window.fa.bg -b - -c 4 -o mean | \
		grep -P "\d$" > ${LOADDIR}/${fraction[${SLURM_ARRAY_TASK_ID}]}'_5mC.window.bg'

	fileName=${fraction[${SLURM_ARRAY_TASK_ID}]}'_5mC.window.bg'
	echo 'Filename is:' $fileName
	sampleName=${fileName%.window.bg}

else 
	echo "Data has been peak called"
	echo 'Changing to directory' $PEAKDIR
	cd ${PEAKDIR}

	## get path to subset directory
	SUBDIR=$(find ${PEAKDIR}/ -name subset)

	## get list of peak called files
	FILES=($(find ${SUBDIR} -name '*.filt'))

	## get name for each datatype and tissue
	sampleName=$(basename ${FILES[${SLURM_ARRAY_TASK_ID}]%.*.filt})

	echo 'Filename is:' ${FILES[${SLURM_ARRAY_TASK_ID}]}

	fileName=$(basename ${FILES[${SLURM_ARRAY_TASK_ID}]}.colsrt)
	echo 'Output sorted file: ' $colsrt

	## create peakfile with only columns 1,2,3,7 in the format:
	# chromosome [s] chrstart [s] chrstop [s] signalValue (measure of total overall enrichment for the region) 
	awk '{print $1,$2,$3,$7}' ${FILES[${SLURM_ARRAY_TASK_ID}]} > ${LOADDIR}/${fileName}
fi



cd ${SCRIPTDIR}
sh ./integrative/segway/processing/genomedataLoad.sh ${sampleName} ${fileName}


## move log files into a folder
cd ${SCRIPTDIR}/integrative/segway/logFiles/${USER}
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}*${SLURM_ARRAY_TASK_ID}* ${SLURM_ARRAY_JOB_ID}