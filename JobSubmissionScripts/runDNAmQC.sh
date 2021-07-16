#!/bin/sh
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=../LogFiles/QCDNAdataMRC.err # error file
#SBATCH --output=../LogFiles/QCDNAdataMRC.log # output file
#SBATCH --job-name=QCDNAdataMRC

## print start date and time
echo Job started on:
date -u

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR

module load Pandoc
module load R/3.6.3-foss-2020a

Rscript performQC.r rmdConfig.mrc 

Rscript -e "rmarkdown::render('QC.rmd', params = list(config='rmdConfig.mrc'), output_file='../../DNAm/QCmetrics/QC.html')"

Rscript -e "rmarkdown::render('QCwithCellType.rmd', params = list(config='rmdConfig.mrc'), output_file='../../DNAm/QCmetrics/QCwithinCellType.html')"


## print finish date and time
echo Job finished on:
date -u