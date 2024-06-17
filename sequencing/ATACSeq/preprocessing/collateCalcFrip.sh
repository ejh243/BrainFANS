#!/bin/bash
## ===================================================================================================================##
##                                ATAC-seq pipeline STEP 3.3: FRiP results                                            ##
## ===================================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/collateCalcFrip.sh <samples>                                      ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script collates the results of peak calling (STEP 3.2) using MACS3 PE at sample level in a       ||
## single csv file.                                                                                                   || 
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <samples> list of samples to collate results for                                                             ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - PEAKDIR/QCOUTPUT/FRIP_all_samples.csv                                                                            ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - bedtools, samtools                                                                                               ||
## - ALIGNEDDIR, PEAKDIR                                                                                              ||
## - .broadPeak.filt for BAMPE.filt.nodup.bam                                                                         ||
##                                                                                                                    ||
## ===================================================================================================================##


SAMPLES=$@

## Output will be a table in a csv file
echo "SampleName,BAMTotalReads,BAMPeaks,ReadsinBAMPeaks" > ${PEAKDIR}/QCOutput/FRIP_all_samples.csv

for sampleName in ${SAMPLES[@]}
do

  ## BAM file used for MACS3 PE peak calling
  bamFile=${sampleName}.filt.nodup.bam 
  
  ## Output sample name
  echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/FRIP_all_samples.csv
  
  if [ -s ${ALIGNEDDIR}/${bamFile} ]
  then 
  	echo -n $(samtools view -c ${ALIGNEDDIR}/${bamFile}), >> ${PEAKDIR}/QCOutput/FRIP_all_samples.csv
  else
  	echo -n NA, >> ${PEAKDIR}/QCOutput/FRIP_all_samples.csv
  fi  
  
  ## MACS3 peaks called from bam files with paired end and number of reads in those peaks
  if [ -s ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt ]
  then
  	echo -n $(wc -l ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/FRIP_all_samples2.csv
  	echo $(bedtools sort -i ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${bamFile} -b stdin -ubam | samtools view -c) >> ${PEAKDIR}/QCOutput/FRIP_all_samples.csv
  else
  	echo -n NA,NA >> ${PEAKDIR}/QCOutput/FRIP_all_samples.csv
  fi

done


