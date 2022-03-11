#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p sq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-193495 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=LogFiles/alignShortRead-%A_%a.o
#SBATCH --error=LogFiles/alignShortRead-%A_%a.e
#SBATCH --job-name=alignShortRead-%A_%a.e
#SBATCH --array=0-3%3 ## runs 50 jobs with 10 at any one time


module load R/3.6.0-foss-2019a


Rscript downloadSynapseData.r ${SLURM_ARRAY_TASK_ID} $1
