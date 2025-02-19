This readme explains how to use the scripts for running the DNAm analysis pipeline

## Prequisites:

1) Copy the provided example config files located in `.../array/DNAm/config` to
a memorable location (ideally next to your data)
	* There is two config files here, one is a `.R` config file and the other is
	`.txt`. The `.txt` config file is to be used when submitting bash scripts,
	whilst the `.R` one is to be used when submitting standalone R scripts.
	* Fill in these config files, providing paths and variable values. Some
	defaults are provided for you
2) Run the setup script found in `array/DNAm/preprocessing/Setup` with bash.
    * Command: `bash .../array/DNAm/preprocessing/Setup path/to/config_file.txt`
		* This config file is the one copied in step 1.
	* With 32GB of RAM and an i7-10610U CPU, the installation process took 1
	hour and 37 minutes. 
	* Compilation of some R packages can halt if your RAM is too low (<1GB),
	consider looking into swapfiles if you think this is happening to you.
3) Create directory structure as outlined in (#directory-structure) where your
data is.
	* Note: This is inline with the defaults provided in the example
	configuration files. If you change the paths for entries such as
	`$METADIR`, your directory structure will need to reflect this change.

## Directory Structure:

```text
Data_directory/
├── 0_metadata/
│   └── sampleSheet.csv # IMPORTANT: This is case sensitive
├── 1_raw/
│   ├── 1234_Grn.idat
│   ├── 1234_Red.idat
│   └── ...
├── 2_gds
└── 3_normalised
```


## OUTPUT:

* html QC report are located in 2_gds/QCmetrics
* text file summarising the QC metrics and filtering are located in 2_gds/QCmetrics
* log files are saved to logFiles in the data directory with the prefix `QCDNAdata`

## Data pre-processing

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
