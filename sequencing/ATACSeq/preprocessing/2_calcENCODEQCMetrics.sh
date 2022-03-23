## calculates ENCODE library complexity metrics (taken from https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#)

## EXECUTION
# sh ./ATACSeq/preprocessing/3_calcENCODEQCMetricsPE.sh <sampleName>
# where 
# <bam file> is the path to a sorted bam file
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR

## REQUIRES the following software
# samtools, samstats, bedtools

## INPUT
# *.filt.nodup.bam aligned, filtered, deduplicated bam file

## OUTPUT
# *.pbc.qc
# *.flagstat.qc

sampleName=$1 

echo "Calculating ENCODE QC metrics"
echo Job started on:
date -u

cd ${ALIGNEDDIR}

mkdir -p ENCODEMetrics

# Run samstats on filtered, dedupped bam file
samtools sort -n --threads 10 ${sampleName}.filt.nodup.bam -O SAM  | SAMstats --sorted_sam_file -  --outf ENCODEMetrics/${sampleName}.flagstat.qc


# =============================
# Compute library complexity
# =============================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics

# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

samtools sort -n ${sampleName}.filt.dupmark.bam -o ${sampleName}.srt.tmp.bam

bedtools bamtobed -bedpe -i ${sampleName}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{if (mt > 0) {printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2} else printf "%d\t%d\t%d\t%d\t%s\t%s\t%s\n",mt,m0,m1,m2,NA,NA,NA }'  > ENCODEMetrics/${sampleName}.pbc.qc

rm ${sampleName}.srt.tmp.bam

if [[ $? == 0 ]]
	then echo "Library complexity metrics calculated:"
	date -u
fi