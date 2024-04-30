#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=BSSeq/logFiles/%u/bismark_conversion-%A.o
#SBATCH --error=BSSeq/logFiles/%u/bismark_conversion-%A.e
#SBATCH --job-name=bismark_conversion-%A.e

## print start date and time
echo Job started on:
date -u

REFGENOME=$1

module load Bismark


bismark_genome_preparation --verbose ${REFGENOME}

echo "Bismark genome prep complete"