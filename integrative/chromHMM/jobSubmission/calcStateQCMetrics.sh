#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=18:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/chromHMM/logFiles/%u/robust-%A.o
#SBATCH --error=integrative/chromHMM/logFiles/%u/robust-%A.e
#SBATCH --job-name=robust

#-----------------------------------------------------------------------#

## print start date and time
echo Job started on:
date -u
	
## needs to be executed from the scripts folder


## load config file provided on command line when submitting job
echo "Loading config file for project: " $1
export INTPROJECT=$1

echo 'learnModel job ID is:' $2
JOBID=$2

source ./integrative/chromHMM/config/config.txt 
echo "Project directory is: " $DATADIR

#-----------------------------------------------------------------------#

cd $MODELDIR
mkdir -p ${QCDIR}

## calculate stability
module purge
module load BEDTools

bedFILES=$(ls *dense.bed )

echo "Generating sorted 200bp split chromosome 1 files"

echo "Files are" $bedFILES

for file in ${bedFILES}
do 
	grep 'chr1	' ${file} | sort -k1,1 -k2,2n | bedtools map -a ${REF}/genome.window.chr1.bg -b - -c 4 | awk '{print $4}'> ${QCDIR}/${file%bed}200.bed
done

cd ${QCDIR}
FILES=$(ls *.200.bed)
paste ${FILES} | column -s $'\t' -t > stateAssignments.bed

#module purge
#module load R/3.6.3-foss-2020a
#
#Rscript #
#
### calculate AIC
#echo "CALCULATE AIC"
#
#cd $LOGDIR
#for i in *.o
#do 
#	j=${i#learnModel-*_} 
#	j=${j%.o}
#	awk '/Writing to file \/lustre\/projects\/Research\_Project\-MRC190311\/integrative\/chromHMM\/[A-Za-z0-9\_]+\/[a-z3-4\_]+\/N\+\_[0-9]+\_segments\.bed/{if (a && a !~ /Writing to file \/lustre\/projects\/Research\_Project\-MRC190311\/integrative\/chromHMM\/[A-Za-z0-9\_]+\/[a-z0-9\_]+\/N\+\_[0-9]+\_segments\.bed/) print a} {a=$2}' $i > ${QCDIR}/likelihood.$j.txt
#done

