This readme explains how to use the scripts for running the DNAm analysis pipeline

PREQUISITES:
	* A config file with the file paths specific to their project 
	* Run the setup script found in .../Setup, provide the config file as a command line argument
	* A file named sampleSheet.csv in the 0_metadata folder which lists the samples to process 
	* idats are in the 1_raw folder

OUTPUT:
	* html QC report are located in 2_gds/QCmetrics
	* text file summarising the QC metrics and filtering are located in 2_gds/QCmetrics
	* log files are saved to logFiles in the data directory with the prefix `QCDNAdata`

#### Data pre-processing

Provided is a SLURM job submission script which automates the preprocessing and can be submitted as follows

`sbatch 1_runDNAmQC.sh /path/to/configFile`
	/path/to/configFile the path to config file which specifics the data and script paths for processing

	* executes loadDataGDS.r ${DATADIR}
	* executes calcQCMetrics.r ${DATADIR} ${REFDIR} [${GENOFILE}]
	* executes clusterCellTypes.r ${DATADIR} ${RCONFIG} 
	* executes Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args ${DATADIR} ${RCONFIG} $USER


During execution log files are temporarily located in the folder you submitted the script from with the filenames QCDNAdata_XXX.log and QCDNAdata_XXX.err, where XXX is replaced with the job number provided by your HPC scheduler. Please ensure you have the correct permissions to create files in this folder, otherwise the job will instantly fail. If the job fails before the script is run to completion, this is where you need to look for troubleshooting. Once the QC pipeline is complete these files are moved to your data directory in a folder called `logFiles`. 


Approach is as follows:

1. QC across all samples
2. QC within cell type
3. Normalisation
4. Characterise dataset
5. Analysis
