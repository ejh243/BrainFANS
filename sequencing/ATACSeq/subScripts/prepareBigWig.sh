#!/bin/bash
## ==========================================================================================================##
##                          PIPPA STEP 2.2: Prepare BigWig files for TSS scores                              ##
## ==========================================================================================================##
## EXECUTION: sbatch ./subScripts/prepareBigWig.sh <sampleName>                                              ||
##                                                                                                           ||
## DESCRIPTION: Produce bigWig file in required format for TSS enrichment score calculation.                 ||
##                                                                                                           ||
## INPUTS:                                                                                                   ||
## $1 -> <sampleName> Name of sample specified in command line. Reads should be previously aligned and       ||
##  duplicates marked                                                                                        ||
##                                                                                                           ||
## OUTPUTS:                                                                                                  ||
## - "sampleName.shifted.bw"                                                                                 ||
##                                                                                                           ||
## REQUIRES:                                                                                                 ||
## - Deeptools, SAMtools (in a conda environment)                                                            ||
## - Variables in config file: ALIGNED_DIR                                                                   ||
## - Duplicated marked, aligned bam file: sampleName.filt.dupmark.bam                                        ||
##                                                                                                           ||
## ==========================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

## sampleName to process
export sampleName=$1
echo "Processing" ${sampleName}

cd ${TMPDIR}

### input BAM (duplicates retained, sorted)
INPUT_BAM="${ALIGNED_DIR}/${sampleName}.filt.dupmark.bam"

# check if bam file exists
if [[ ! -f "${INPUT_BAM}" ]]
then 
  { echo "Aligned, duplicated marked, reads file for ${sampleName} not found. Please run STEP 1.3 ALIGN first." ; exit 1; }
else
  echo "Found files: ${INPUT_BAM}"
fi

## =========== ##
##  BAM to BW  ##
## =========== ##

echo "Sorting sample"
SORTED_DUP_BAM="${ALIGNED_DIR}/${sampleName}_filt.dup.sorted.bam"

### Sort bam file if not sorted
samtools sort "${INPUT_BAM}" > "${SORTED_DUP_BAM}"
samtools index "${SORTED_DUP_BAM}"

### ---- Tn5 shift using bamCoverage offsets ----
# + strand (+4 bp)
bamCoverage -b "${SORTED_DUP_BAM}" -o "${TSS_DIR}/${sampleName}.plus.bw" \
  --filterRNAstrand forward \
  --Offset 4 \
  --normalizeUsing CPM \
  --binSize 10 \
  --numberOfProcessors 1

# - strand (-5 bp)
bamCoverage -b "${SORTED_DUP_BAM}" -o "${TSS_DIR}/${sampleName}.minus.bw" \
  --filterRNAstrand reverse \
  --Offset -5 \
  --normalizeUsing CPM \
  --binSize 10 \
  --numberOfProcessors 1

### ---- merge strands ----
bigwigCompare -b1 "${TSS_DIR}/${sampleName}.plus.bw" -b2 "${TSS_DIR}/${sampleName}.minus.bw" \
  --operation add \
  -o "${TSS_DIR}/${sampleName}.shifted.bw"

## remove intermediate files
rm -f "${TSS_DIR}/${sampleName}.plus.bw" "${TSS_DIR}/${sampleName}.minus.bw" "${SORTED_DUP_BAM}"

if [[ ! -f "${TSS_DIR}/${sampleName}.shifted.bw" ]]
then 
  { echo "Shifted bigWig file for ${sampleName} could not be completed." ; exit 1; }
else
  echo "Shifted bigWig file for ${sampleName} is done."
  echo Job finished on:
  date -u 
fi 