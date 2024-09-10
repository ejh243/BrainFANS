#!/bin/bash
## calculates ENCODE library complexity metrics (taken from https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#)

## EXECUTION
# sh ./BSSeq/preprocessing/2_calcENCODEQCMetrics.sh <sampleName>
# where 
# <bam file> is the path to a sorted bam file
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR

## REQUIRES the following software
# samtools, samstats, bedtools

## INPUT
# *.bam aligned, filtered, deduplicated bam file
# *.filt.dupmark.bam filtered, duplicate-marked bam file

## OUTPUT
# *.pbc.qc
# *.flagstat.qc

sampleName=$1 
echo
echo "Calculating ENCODE QC metrics"
echo Job started on:
date -u

cd ${ALIGNEDDIR}

## compute average coverage across genome
## each replicate should have 30X coverage.

echo 'Computing average coverage'
samtools sort ${sampleName}*nodup.bam | bedtools genomecov -ibam stdin | \
	awk '{SUM+=$2 * $3}{ a[$4]++ } END { for (b in a) { TOT+=b } {print SUM/TOT} }' > ENCODEMetrics/${sampleName}.qc


## conversion efficiency
# the C to T conversion rate should be ≥98%

echo "Computing conversion efficiency"
if test -d spikeAlignments 
then 
	var=$( grep 'C methylated in CpG context:' spikeAlignments/${sampleName}*.txt | awk '{ print $6 }' | rev | cut -c2- | rev )
	echo "100 - $var" | bc >> ENCODEMetrics/${sampleName}.qc
else
	echo 'No separate alignment found. Conversion efficiency not calculated'
fi


## Compute correlation part 1

#cd $METHYLDIR

#zcat ${sampleName}*bismark.cov.gz | grep "chr1	" > ENCODEMetrics/${sampleName}.chr1.cov.bg

if [[ $? == 0 ]]
	then echo "Part 1 ENCODE metrics calculated"
	date -u
fi
