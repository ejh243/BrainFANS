#!/bin/bash
## calculates number of reads within peaks for sample level peaks

## EXECUTION
# sh ./ChIPSeq/preprocessing/6_calcFrip.sh <sampleName>
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# bedtools, samtools

## INPUT
# blacklist filtered peak lists
#

## OUTPUT
# ${PEAKDIR}/QCOutput/FRIP_*.csv

sampleName=$1
echo
echo "Starting calculating frip for" ${sampleName} "at: "
date -u


## save output in txt file
echo "SampleName,BAMTotalReads,MACS2PEPeaks,ReadsinMACS2PEPeaks" > ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv

echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv

# formulate peak input filenames

# count total reads in peak input files
if [ -s ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam ]
then 
	echo -n $(samtools view -c ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam), >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
else
	echo -n NA, >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
fi

# count reads in peak regions

## MACS2 peaks called from bam files with paired end flag

if [ -s "${PEAKDIR}/${sampleName}.narrowPeak.filt" ] || [ -s "${PEAKDIR}/${sampleName}.broadPeak.filt"  ]
then
	echo -n $(wc -l ${PEAKDIR}/${sampleName}*Peak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
	echo $(bedtools sort -i ${PEAKDIR}/${sampleName}*Peak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam -b stdin -ubam | samtools view -c) >> ${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv
else
	echo NA,NA >> "${PEAKDIR}/QCOutput/FRIP_${sampleName}.csv"
fi
