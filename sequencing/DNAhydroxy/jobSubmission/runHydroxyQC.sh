#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/CEGX5hmCQC.err # error file
#SBATCH --output=LogFiles/CEGX5hmCQC.log # output file
#SBATCH --job-name=CEGX5hmCQC



## print start date and time
echo Job started on:
date -u

## needs to be executed from the scripts folder

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

## load config parameters
source hydroxy/CGEX/config.txt

echo Starting peak calling at:
date -u

module load MACS2/2.1.2.1-foss-2017b-Python-2.7.14

cd hydroxy/CGEX/

## For QC purposes first run peak calling on each sample
./macsPeakCalling.sh

## Calc QC metrics

## need indexed bamfiles
BAMFILES=($(ls ${ALIGNEDDIR}/*_L00.bml.GRCh38.karyo.deduplicated.bam))
module load SAMtools
for f in ${BAMFILES[@]}
do
	samtools index ${f}
done

## Run QC report
module load Pandoc/2.5
module load R/3.6.3-foss-2020a

Rscript calcQCMetrics.r config.r 

Rscript -e "rmarkdown::render('QCReport.rmd', params = list(config='config.r'), output_file='QCJuly2021.html')"


## print finish date and time
echo Job finished on:
date -u