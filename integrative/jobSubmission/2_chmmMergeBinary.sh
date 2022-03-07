#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=25:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=5 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=integrative/logFiles/mergeBinary-%A.o
#SBATCH --error=integrative/logFiles/mergeBinary-%A.e
#SBATCH --job-name=mergeBinary


# This performs mergebinary on the 
# requires all data-type/project-name directories to be merged on the command line. i.e.
# sbatch /path/to/jobscript/ ATACSeq/rizzardi WGBS/rizzardi ChIPSeq/epiGaba DNAhydroxySeq/mrcScz 

## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder
echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR
echo 'Loading config file'
source ./integrative/config/config.txt

##shift # Cause all of the positional parameters $2 to $n to be renamed $1 to $(n-1), and $1 to be lost.


DATADIRS=($@)

TYPES=( "${DATADIRS[@]%/*}" )
echo ${TYPES[@]}

#create full paths to binarised data from command line arguements
DATADIRS=( "${DATADIRS[@]/#/${PROJDIR}/}" )  # prepend the path to main project directory folder to the data-type/project variable
DATADIRS=( "${DATADIRS[@]/%//3_aligned/binarisedData}" ) # append the path to the binarisation



mkdir -p ${MERGEDIR}/1_input ${MERGEDIR}/2_output

##create softlinks of all data to a merged folder
for (( x=0; x<${#DATADIRS[@]}; x++ ))
do
	if [ -d ${MERGEDIR}/input/${TYPES[${x}]} ] #check if data type directory exists
	then 
		echo "Softlinks already created for data type: " ${TYPES[${x}]}
	else
		echo "Creating link to " ${DATADIRS[${x}]}
		echo
		ln -s ${DATADIRS[${x}]} ${MERGEDIR}/input/${TYPES[${x}]}
	fi
done

echo "Running MergeBinary"

module load Java

java -mx2400M -jar ${CHROMHMM}/ChromHMM.jar MergeBinary ${MERGEDIR}/input ${MERGEDIR}/2_output

if [[ $? == 1 ]]
then
	echo "Finished merge"
fi