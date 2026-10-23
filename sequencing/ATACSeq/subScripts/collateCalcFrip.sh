#!/bin/bash
## ===================================================================================================================##
##                                           PIPPA STEP 3.3: FRiP results                                             ##
## ===================================================================================================================##
## EXECUTION: sbatch ./subScripts/collateCalcFrip.sh                                                                  ||
##                                                                                                                    ||
## DESCRIPTION: This script collates the results of peak calling (STEP 3.2 PEAKS) using MACS3 PE at sample level in a ||
## single csv file.                                                                                                   || 
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - PEAK_DIR/QCOUTPUT/FRIP_all_samples.csv                                                                           ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - BEDTools, SAMtools (in a conda environment)                                                                      ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR,PEAK_DIR_SAMPLES                                                 ||
## - Peak files for samples: <sampleName>.narrowPeak.filt                                                             ||
## - Bam files for samples: <sampleName>.filt.nodup.bam                                                               ||
##                                                                                                                    ||
## ===================================================================================================================##


## ============ ##
##    SET UP    ##
## ============ ##

mapfile -t SAMPLES < "${SAMPLES_LIST}"
echo "Total number of samples: "" ""${#SAMPLES[@]}"""

FRIP_FILE="${PEAK_DIR}/QCOutput/FRIP_all_samples.csv"

header="SampleName,BAMTotalReads,BAMPeaks,ReadsinBAMPeaks"

## ============ ##
##   COLLATE    ##
## ============ ##

if [[ ! -s "${FRIP_FILE}" ]]; then
    printf "%b\n" "${header}" > "${FRIP_FILE}"
fi

{ 
  for sampleName in "${SAMPLES[@]}"
  do
    if grep -q "^${sampleName}," "${FRIP_FILE}"; then
      continue
    fi
  
    ## BAM file used for MACS3 PE peak calling
    INPUT_BAM="${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"
    
    ## Output sample name
    echo -n ${sampleName},
    
    if [ -s "${INPUT_BAM}" ]
    then 
    	echo -n $(samtools view -c "${INPUT_BAM}"), 
    else
    	echo -n NA, 
    fi  
    
    ## MACS3 peaks called from bam files with paired end and number of reads in those peaks
    INPUT_PEAKS="${PEAK_DIR_SAMPLES}/${sampleName}.narrowPeak.filt"
    
    if [ -s "${INPUT_PEAKS}" ]
    then
    	echo -n $(wc -l "${INPUT_PEAKS}" | cut -f1 -d' '), 
    	echo $(bedtools sort -i "${INPUT_PEAKS}" | bedtools merge -i stdin | bedtools intersect -u -a "${INPUT_BAM}" -b stdin -ubam | samtools view -c) 
    else
    	echo -n NA, 
      echo NA 
    fi
  done

} >> "${FRIP_FILE}"


if [[ ! -f "${FRIP_FILE}" ]]
then
  { echo "Peak results could not be collated. Please check STEP 3.2 PEAKS." ; exit 1;}
else
  echo "Peak calling results have been collated."
  echo Job finished on:
  date -u
fi
