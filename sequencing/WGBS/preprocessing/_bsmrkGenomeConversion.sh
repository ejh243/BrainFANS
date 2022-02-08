#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=48:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=LogFiles/WGBS/bismark_conversion-%A.o
#SBATCH --error=LogFiles/WGBS/bismark_conversion-%A.e
#SBATCH --job-name=bismark_conversion-%A.e

## print start date and time
echo Job started on:
date -u

module load Bismark

bismark_genome_preparation --verbose /gpfs/mrc0/projects/Research_Project-MRC190311/References/grch38/ncbi/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta

echo "Bismark genome prep complete"