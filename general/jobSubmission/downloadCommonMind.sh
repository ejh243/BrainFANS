#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --output=general/logFiles/downloadCommonmind-%A_%a.o
#SBATCH --error=general/logFiles/downloadCommonmind-%A_%a.e
#SBATCH --job-name=downloadCommonmind-%A_%a.e
#SBATCH --array=5-10 


module load R/3.6.0-foss-2019a


Rscript general/processing/downloadCommonMind.r ${SLURM_ARRAY_TASK_ID} $1
