#!/bin/bash 
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=1:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=1 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=general/logFiles/filePermission-%A.o
#SBATCH --error=general/logFiles/filePermission-%A.o
#SBATCH --job-name=filePermission

## to submit sbatch /path/to/jobscript /datatypedir/to/change/within

## print start date and time
echo Job started on:
date -u

# repository directory
REPO_DIR=$1

echo "ATTENTION! Will change all files bar those in 1_raw to read/write/execute permissions for those in your group"
echo
echo "--- Setting rwx permissions to group on all non-raw data files in ${REPO_DIR} owned by ${USER}"


if [ $# == 1 ];
then
    chmod 774 $(find ${REPO_DIR} -path **1_raw -prune -o -user ${USER} -print)
elif [[ $2 =~ 'END' ]];
then
    echo "ATTENTION! This will also change all files in 1_raw to read/write/execute permissions for those in your group"
    files=($(find ${REPO_DIR} -path -prune -o -user ${USER} -print))
    for file in ${files[@]}
    chmod 774 $file
fi




echo "--- Done"
echo "---";