#!/bin/bash
## ==========================================================================================================##
##                              PIPPA STEP 2.3: Calculate TSS enrichment score                               ##
## ==========================================================================================================##
## EXECUTION: sbatch ./subScripts/calcTSSEnrich.sh <sampleName>                                              ||
##                                                                                                           ||
## DESCRIPTION: This calculates the TSS enrichment score of a selected sample.                               ||
##                                                                                                           ||
## INPUTS:                                                                                                   ||
## $1 -> <sampleName> Name of sample specified in command line.                                              ||
##                                                                                                           ||
## OUTPUTS:                                                                                                  ||
## - ALIGNED_DIR/bigWig/tssEnrich/${sampleName}_TSS.mat.gz                                                   ||
## - ALIGNED_DIR/bigWig/tssEnrich/${sampleName}.profile.txt                                                  ||
## - ALIGNED_DIR/bigWig/tssEnrich/${sampleName}.profile.pdf                                                  ||
##                                                                                                           ||
## REQUIRES:                                                                                                 ||
## - deeptools (in a conda environment)                                                                      ||
## - Variables in config file: ALIGNED_DIR, REFTSS                                                           ||
## - Shifted BW file for sampleName: <sampleName>.shifted.bw                                                 ||
## - Reference TSS file (path specified in config.txt file)                                                  ||
##                                                                                                           ||
## ==========================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## sampleName to process
export sampleName=$1
echo "Processing" ${sampleName}
cd ${TMPDIR}

INPUT_BW="${TSS_DIR}/${sampleName}.shifted.bw"

# check if bigWig file exists
if [[ ! -f "${INPUT_BW}" ]]
then 
  { echo "BigWig file for ${sampleName} not found. Please run STEP 2.2 BW first." ; exit 1; }
else
  echo "Found files: ${INPUT_BAM}"
fi

# check if reference TSS file has been provided
if [[ ! -f "${REFTSS}" ]]
then 
  { echo "TSS reference file not found. Please specify this in the config.txt file." ; exit 1; }
else
  echo "${REFTSS} is reference file used"
fi

## ===================== ##
##  CALCULATE TSS score  ##
## ===================== ##

OUTPUT_MAT="${TSS_DIR}/${sampleName}_TSS.mat.gz"
 
if [[ ! -f "${OUTPUT_MAT}" ]]; then
  computeMatrix reference-point \
    --referencePoint TSS \
    -b 2000 -a 2000 \
    -R "${REFTSS}" \
    -S "${INPUT_BW}" \
    --skipZeros -o "${OUTPUT_MAT}"
fi 

if [[ ! -f "${OUTPUT_MAT}" ]]
then 
  { echo "File with TSS score ${OUTPUT_MAT} not created. Check STEP 2.2" ; exit 1; }
else
  echo "TSS score for ${sampleName} is calculated."
  echo Job finished on:
  date -u 
fi 