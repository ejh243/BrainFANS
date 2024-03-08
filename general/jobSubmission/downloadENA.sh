#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=100:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --mem=50G
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=general/logFiles/downloadPsych-%A.o
#SBATCH --error=general/logFiles/downloadPsych-%A.e
#SBATCH --job-name=downloadPsych


## if raw data directory is empty, download files from ENA with specified ftp list
echo "Downloading files"
ACCLIST=${METADIR}/file*
python ${SCRIPTDIR}/general/processing/ftp_url.py $ACCLIST #generate ena download list
cd ${SCRIPTDIR}
sh ./general/jobSubmission/downloadFromENA.sh  
fi

wget -i ${DATADIRPE}/ftp_url_list.txt