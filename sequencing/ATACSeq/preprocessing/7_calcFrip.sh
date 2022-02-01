## calculates number of reads within peaks for sample level peaks

## EXECUTION
# sh ./ATACSeq/preprocessing/7_calcFrip.sh <sampleName>
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
# ${PEAKDIR}/QCoutput/FRIP_*.csv


sampleName=$1

## save output in txt file
echo "SampleName,TagAlignTotalReads,BAMTotalReads,MACS2TagAlignPeaks,ReadsinMACS2TagAlignPeaks,MACS2PEPeaks,ReadsinMACS2PEPeaks" > ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv

echo -n ${sampleName}, >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv

# formulate peak input filenames
tagAlign=${sampleName}.tn5.tagAlign
bam=${sampleName}.filt.nmsrt.nodup.bam # NB identical to ${sampleName}_depDup_q30.bam

# count total reads in peak input files
if [ -s ${ALIGNEDDIR}/${tagAlign} ]
then 
	echo -n $(zcat ${ALIGNEDDIR}/${tagAlign} | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
else 
	echo -n NA, >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
fi
if [ -s ${ALIGNEDDIR}/${bam} ]
then 
	echo -n $(samtools view -c ${ALIGNEDDIR}/${bam}), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
else
	echo -n NA, >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
fi

# count reads in peak regions

## MACS2 peaks called from shifted tag align file
if [ -s ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak.filt ]
then
	echo -n $(wc -l ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
	echo -n $(bedtools sort -i ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${tagAlign} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
else
	echo -n NA,NA, >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
fi

## MACS2 peaks called from bam files with paired end flag
if [ -s ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak ]
then
	echo -n $(wc -l ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
	echo $(bedtools sort -i ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${bam} -b stdin -ubam | samtools view -c) >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
else
	echo NA,NA >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
fi
