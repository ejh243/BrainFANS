#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=96:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=integrative/chromTools/logFiles/%u/complete-%A_%a.e # error file
#SBATCH --output=integrative/chromTools/logFiles/%u/complete-%A_%a.o # output file
#SBATCH --job-name=complete

#-----------------------------------------------------------------------#

## needs to be run with array number 0-number of targets 
## at the moment can only be run from $USER=jms260 !!! (hard coded chromtools path)

## print start date and time
echo Job started on:
date -u

PROJECT=$1
source ./integrative/chromTools/config/config.txt 

## check directories
echo "Project directory is: " $DATADIR
echo 'Script directory is: ' ${SCRIPTDIR}

##check array specified and exit if not
if [[ ${SLURM_ARRAY_TASK_ID} == '' ]]
then 
    { echo "Job does not appear to be an array. Please specify --array on the command line." ; exit 1; }
fi

#-----------------------------------------------------------------------#

## column of targets
targetCol=$(awk '$1 == "target"{print NR;exit} ' RS="," ${METADIR}/sampleSheet.csv)
sampleCol=$(awk '$1 == "sampleID"{print NR;exit} ' RS="," ${METADIR}/sampleSheet.csv)
controlCol=$(awk '$1 == "controlID"{print NR;exit} ' RS="," ${METADIR}/sampleSheet.csv)

## get targets
IFS=$'\n' # need to set this as \n rather than default - space, \t and then \n - so that elements are expanded using \n as delimiter
target=($(awk -F"," -v col="$targetCol" '(NR>1) {print $col}' ${METADIR}/sampleSheet.csv | sort -u))
unset IFS

target=${target[${SLURM_ARRAY_TASK_ID}]}
echo "Target is: " $target


SAMPLEIDS=($(grep "${target}" ${METADIR}/sampleSheet.csv | awk -F',' -v col="$sampleCol" '{print $col}'))

cd $ALIGNEDDIR

module load BEDTools

## convert to BED file

echo "Sample IDs are" ${SAMPLEIDS[@]}

if [ ! -z "$controlCol" ]; #if there are controls (if controlCol is not empty)
then
    CONTROLIDS=($(grep "${target}" ${METADIR}/sampleSheet.csv | awk -F',' -v col="$controlCol" '{print $col}'))

    for sampleName in ${CONTROLIDS[@]};
    do
        if [ ! -f ${sampleName}.bed ];
        then
            echo "converting" $sampleName "to BED"
            bamToBed -i ${sampleName}.filt.nodup.bam > ${sampleName}.bed
            echo $sampleName
        fi
    done

    ctrlfiles=${CONTROLIDS[@]/%/.bed}
fi


for sampleName in ${SAMPLEIDS[@]};
do
    if [ ! -f ${sampleName}.bed ];
    then
        bamToBed -i ${sampleName}.filt.nodup.bam > ${sampleName}.bed
        echo $sampleName
    fi
done

module purge
module load Python/3.7.2-GCCcore-6.4.0
source ~/.chromTools/bin/activate

files=${SAMPLEIDS[@]/%/.bed}

echo "Files are: " $files

CHROMTOOLDIR=${CHROMTOOLDIR}

mkdir -p $CHROMTOOLDIR


if [ ! -z "$controlCol" ];
then
    chromTools complete -f $files -c $ctrlfiles -o ${CHROMTOOLDIR}/${target} --force-overwrite
else
    chromTools complete -f $files -o ${CHROMTOOLDIR}/${target} --force-overwrite
fi



## move log files into a folder
cd ${SCRIPTDIR}/integrative/chromTools/logFiles/$USER
mkdir -p ${SLURM_ARRAY_JOB_ID}
mv *${SLURM_ARRAY_JOB_ID}* ${SLURM_ARRAY_JOB_ID}