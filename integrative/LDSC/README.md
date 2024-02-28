# LDSC scripts

The scripts in this folder are desinged to create an annotation file from peak called data (from ATAC, ChIP *etc.*), then uses [ldsc](https://github.com/bulik/ldsc) to find the partitioned heritability for the dataset (partitioned by said peaks).

## File structure
The basic file structure that is expected with these scripts is as follows:
```text
Main directory for ldsc files
├── Reference files
│   ├── weights
│   ├── PLINK files
│   └── frq files
├── Annotation files
├── SNP lists
├── Collection of GWAS traits
└── Outputs
    ├── Reference ldsc files from PLINK files (under annotations)
    └── Heritability analysis outputs

Peak calls data directory
├── 0_metadata
├── 1_raw
├── 2_trimmed
├── 3_aligned
├── 4-called-peaks
│   └── organised peak files
└── [Downstream analysis]
```

## Configuration file

Please save all configuration files in the same directory, ideally somewhere near the ldsc files (though this is not a requirement).
These scripts require the following configuration files:

### config.txt
```text
# File paths in data directory
export MAIN_DIR=path/to/main/data/directory
export LD_REFERENCE_DIR=${MAIN_DIR}/path/to/reference/files
export LD_ANNOTATION_DIR=${MAIN_DIR}/path/to/annot/files
export LD_GWAS_TRAITS_DIR=${MAIN_DIR}/path/to/gwas/traits
export SNP_LISTS_DIR=${MAIN_DIR}/path/to/snp/lists
export OUTPUTS_DIR=${MAIN_DIR}/path/to/outputs

# (optional) Pattern for GWAS traits to consider, default is *
# The file extensions are added for you, do not include these in this variable
export GWAS_PATTERN=[glob pattern]

# Directory to put log files from scripts into
export LOG_DIR=path/to/log/files

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

# LD annotation directory (make sure this is the same as the directory pointed
# to in config.txt)
ld_annotation_dir <- "path/to/ld/annotation/files"
```