This readme explains how to use the scripts for running the DNAm analysis pipeline

PREQUISITES:
	* Scripts submitted from the scripts/array folder
	* Unix config file located scripts/array/DNAm/config/config.txt that generates the file paths specific to the project 
	* Users require a folder named with their username in the DNAm/logFiles directory 
	* A file named sampleSheet.csv in the 0_metadata folder which lists the samples to process 
	* idats are in the 1_raw folder

OUTPUT:
	* html QC reports are located in 2_gds/QCmetrics
	* text file summarising the QC metrics and filtering are located in 2_gds/QCmetrics

#### Data pre-processing:

Parameters in [] are optional.

1. sbatch DNAm/jobSubmission/1_runDNAmQC.sh <project-name> <mkdownConfig>
	<project-name> the name of the project folder in the DNAm data folder
	<mkdownConfig> path to txt file with parameters for normalisation and which qc checks to perform. Note needs to be absolute path.

	* executes DNAm/preprocessing/loadDataGDS.r ${DATADIR}
	* executes DNAm/preprocessing/calcQCMetrics.r ${DATADIR} ${REFDIR} [${GENOFILE}]
	* executes Rscript -e "rmarkdown::render('DNAm/preprocessing/QC.rmd', output_file='QC.html')" --args ${DATADIR} <mkdownConfig> $USER
	* executesDNAm/preprocessing/clusterCellTypes.r ${DATADIR} <mkdownConfig> 
	* Rscript -e "rmarkdown::render('DNAm/preprocessing/QCwithinCellType.rmd', output_file='QCwithinCellType.html')" --args ${DATADIR} <mkdownConfig> $USER
	



Approach is as follows:

1. QC across all samples
2. QC within cell type
3. Normalisation
4. Characterise dataset
5. Analysis