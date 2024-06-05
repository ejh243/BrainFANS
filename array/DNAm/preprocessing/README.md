This readme explains how to use the scripts for running the DNAm analysis pipeline

PREQUISITES:
	* A config file with the file paths specific to their project 
	* Users require a folder named with their username in the DNAm/logFiles directory 
	* A file named sampleSheet.csv in the 0_metadata folder which lists the samples to process 
	* idats are in the 1_raw folder

OUTPUT:
	* html QC reports are located in 2_gds/QCmetrics
	* text file summarising the QC metrics and filtering are located in 2_gds/QCmetrics

#### Data pre-processing

Provided is a SLURM job submission script which automates the preprocessing and can be submitted as follows

`sbatch 1_runDNAmQC.sh /path/to/configFile`
	/path/to/configFile the path to config file which specifics the data and script paths for processing

	* executes loadDataGDS.r ${DATADIR}
	* executes calcQCMetrics.r ${DATADIR} ${REFDIR} [${GENOFILE}]
	* executes Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args ${DATADIR} ${RCONFIG} $USER
	* executes DNAm/preprocessing/clusterCellTypes.r ${DATADIR} ${RCONFIG} 
	* Rscript -e "rmarkdown::render('QCwithinCellType.rmd', output_file='QCwithinCellType.html')" --args ${DATADIR} ${RCONFIG} $USER
	



Approach is as follows:

1. QC across all samples
2. QC within cell type
3. Normalisation
4. Characterise dataset
5. Analysis
