# How to use the DNAm preprocessing(QC) pipeline

## Required Software

The DNAm QC pipeline requires the following software in order to run. Ensure
that these software can be found on your `PATH` environment variable, *i.e.*
you can run these from the command line directly. You can check this by using
`which [command]`, for example `which conda`.

- [bash](https://www.gnu.org/software/bash/) (>=4.2.46(2))
- [SLURM Workload Manager](https://slurm.schedmd.com/overview.html) (>=20.02.3)
- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) (>=v23.10.0)

All additional dependencies are downloaded and installed via conda upon
completing [setup](#setup).

## Setup

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
3) Create the expected [directory structure](#directory-structure) where your
data is.
	* Note: This is inline with the defaults provided in the example
	configuration files. If you change the paths for entries such as
	`$METADIR`, your directory structure will need to reflect this change.

## Directory Structure

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

## Running the preprocessing(QC) pipeline

A bash script is provided in `.../array/DNAm/preprocessing/jobSubmission` which
runs the full pipeline for you. It is to be submitted using SLURM Workload
Manager with the `.txt` config file provided as the first (and only) argument.

Example:

```bash
sbatch 1_runDNAmQC.sh /path/to/configFile
```

This script executes the following R scripts for you, providing the necessary
command line arguments using the provided config file:

* `checkConfigRFile.r`
* `checkColnamesSampleSheet.r`
* `loadDataGDS.r` 
* `calcQCMetrics.r` 
* `clusterCellTypes.r`
* Rendering of `QC.rmd`
* `normalisation.r`
* `CETYGOdeconvolution.r`

## Outputs

* A [gds file](https://pmc.ncbi.nlm.nih.gov/articles/PMC6419913/) (and a
normalised gds) file containing summarised information from your input IDAT
files.
* A html QC report located in 2_gds/QCmetrics
* Comma separated value (csv) files summarising the QC metrics and filtering
are located in 2_gds/QCmetrics
* Rdata files containing information on normalised data and QC metrics

## Log files

All log files are saved to `logFiles` in the data directory. These log files
will have a name of the form `QCDNAdata_%j.err` and `QCDNAdata_%j.log`. Where:

* `%j` is the job ID given to the job by SLURM Workload Manager
* `.err` contains any text that was sent to the standard error stream
	* This doesn't necessarily only contain errors. Warning messages and a lot
	of R output is written here (wrongfully, in my opinion)
* `.log` contains any text that was sent to the standard output stream

If the script fails in some way (job cancelled by user or a crash occurs), the
log files will be located in the directory you submitted the script from. This
is a quirk that is hard to avoid with job schedulers, if you can't find log
files in the expected `logFiles` directory, this is likely why.
