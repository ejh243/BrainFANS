#!/bin/bash
## ===================================================================================================================##
##                                ATAC-seq pipeline STEP 3.3: FRiP results                                            ##
## ===================================================================================================================##
## EXECUTION: sh ./subScripts/collateCalcFrip.sh <array-samples>                                                      ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script collates the results of peak calling (STEP 3.2 PEAKS) using MACS3 PE at sample level in a ||
## single csv file.                                                                                                   || 
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <array-samples> array of samples to collate results for                                                      ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - PEAK_DIR/QCOUTPUT/FRIP_all_samples.csv                                                                           ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - BEDTools, SAMtools                                                                                               ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR                                                                  ||
## - Peak files for samples: <sampleName>.broadPeak.filt                                                              ||
## - Bam files for samples: <sampleName>.filt.nodup.bam                                                               ||
##                                                                                                                    ||
## ===================================================================================================================##

echo Job started on:
date -u
  
SAMPLES=$@
mkdir -p ${PEAK_DIR}/QCOutput

## Output will be a table in a csv file
echo "SampleName,BAMTotalReads,BAMPeaks,ReadsinBAMPeaks" > ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv

for sampleName in ${SAMPLES[@]}
do

  ## BAM file used for MACS3 PE peak calling
  bamFile=${sampleName}.filt.nodup.bam 
  
  ## Output sample name
  echo -n ${sampleName}, >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
  
  if [ -s ${ALIGNED_DIR}/${bamFile} ]
  then 
  	echo -n $(samtools view -c ${ALIGNED_DIR}/${bamFile}), >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
  else
  	echo -n NA, >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
  fi  
  
  ## MACS3 peaks called from bam files with paired end and number of reads in those peaks
  if [ -s ${PEAK_DIR}/BAMPE/${sampleName}.broadPeak.filt ]
  then
  	echo -n $(wc -l ${PEAK_DIR}/BAMPE/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
  	echo $(bedtools sort -i ${PEAK_DIR}/BAMPE/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNED_DIR}/${bamFile} -b stdin -ubam | samtools view -c) >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
  else
  	echo -n NA, >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
    echo NA >> ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv
  fi
done

if [[ ! -f ${PEAK_DIR}/QCOutput/FRIP_all_samples.csv ]]
then
  { echo "Peak results could not be collated. Please check STEP 3.2 PEAKS." ; exit 1;}
else
  echo "Peak calling results have been collated."
  echo Job finished on:
  date -u
fi

