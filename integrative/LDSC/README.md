# LDSC scripts

The scripts in this folder are desinged to create an annotation file from peak called data (from ATAC, ChIP *etc.*), then uses [ldsc](https://github.com/bulik/ldsc) to find the partitioned heritability for the dataset (partitioned by said peaks).

In order to use the ldsc tool, you need to create a conda environment as per the official [documentation](https://github.com/bulik/ldsc?tab=readme-ov-file#getting-started). Please take note of the location of this conda environment so that it can be put into the configuration file. To find the location of your conda environments, use `conda env list` in the terminal.

GWAS traits are expected to be in sumstats.gz format. If this is not the case, please use ldsc to convert any GWAS analysis text files into this file format (munge_stats.py). The file format for sumstats is given [here](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format).

## File structure
The basic file structure that is expected with these scripts is as follows:
```text
Main directory for ldsc files
├── reference_files
│   ├── weights (.l2.ldscore.gz)
│   ├── plink_files (.bim, .fam, .bed)
│   └── frq_files
├── annotation_files
│   └── all annotation files (.annot.gz)
├── SNP lists (.snp)
├── gwas_traits (.sumstats.gz)
└── Outputs
    ├── annotation_prefix_1
    │    ├── ld scores under annotation file used (.l2.ldscore.gz, .log etc.)
    │    └── Heritability analysis outputs (.results, .log)
    └── annotation_prefix_2
         ├── ld scores under annotation file used (.l2.ldscore.gz, .log etc.)
         └── Heritability analysis outputs (.results, .log)


Peak calls data directory
├── 0_metadata
├── 1_raw
├── 2_trimmed
├── 3_aligned
├── 4-called-peaks
│   └── organised peak files (to be used in annotation file creation)
└── [Downstream analysis]
```

## Configuration file

Please save all configuration files in the same directory, ideally somewhere near the data files for the project (though this is not a requirement).
You will need to parse the full (or relative) path to the folder that contains these configuration files when calling runLDSCNeuralPeaks.sh

### config.txt
```text
# File paths in data directory
export MAIN_DIR=path/to/main/data/directory
export LD_REFERENCE_DIR=${MAIN_DIR}/reference_files
export LD_ANNOTATION_DIR=${MAIN_DIR}/annotation_files
export LD_GWAS_TRAITS_DIR=${MAIN_DIR}/gwas_traits
export SNP_LISTS_DIR=${MAIN_DIR}/path/to/snp/lists
export OUTPUTS_DIR=${MAIN_DIR}/path/to/outputs

# File prefixes in data directory, prefix of the file up until the component
# before chromosome. Example: 1000G.EUR.QC
# Note that your frq and plink files are assumed to have the same prefix
export REFERENCE_PREFIX="prefix"
export ANNOTATION_PREFIX="prefix"
export SNP_LIST_PREFIX="prefix"
export WEIGHTS_PREFIX="prefix"

# (optional) Pattern for GWAS traits to consider, default is *
# The file extensions are added for you, do not include these in this variable
export GWAS_PATTERN=[glob pattern]

# Directory to put log files from scripts into
export LOG_DIR=path/to/log/files

# Conda functionality
export LDSC_CONDA_ENVIRONMENT=path/to/conda/environment
export CONDA_SHELL=path/to/conda/etc/directory # found in conda installation

# ldsc software file path
export LD_SOFTWARE_DIR=path/to/ldsc/root/directory
```

### config.R

```R
# paths to projects considered
atac_project_path <- "path/to/atac/project"
cegx_project_path <- "path/to/cegx/project"
# etc.

# Vector of all peak files to call in each project
peaks_file_paths <-
  c(paste0(atac_project_path, "path/to/neun/gapped/peak"),
    paste0(atac_project_path, "path/to/sox10/gapped/peak"),
    paste0(cegx_project_path, "path/to/narrow/peak"),
    ...)

# Corresponding names for each of the peak files
peaks_names <-
  c("NeuN_gapped_peak", "SOX10_gapped_peak", "CEGX_narrow_peak", ...)

# Input files, bims (PLINK files) and Baseline annotations.
# Specify the full file path and then the prefix of the file names (up to the 
# chromosome number, e.g. 1000G.EUR.hg38.)
# Don't forget to include the "." before the chromosome number.
baseline_file_prefix <- "path/to/baseline/directory/baseline_file_prefix"
bim_file_prefix <- "path/to/plink/files/directory/plink_file_prefix"

# LD annotation directory (make sure this is the same as the directory pointed
# to in config.txt)
ld_annotation_dir <- "path/to/ld/annotation/files"

# This is the prefix to use with the annotation files.
ld_annotation_prefix <- "prefix_name"
```