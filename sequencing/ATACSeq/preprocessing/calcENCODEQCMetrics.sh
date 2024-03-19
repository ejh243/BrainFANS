## ===================================================================================================================##
##                    ATAC-seq pipeline STEP 1.4: Pre-analysis -- ENCODEQC metrics calculation                        ##
## ===================================================================================================================##
## EXECUTION: sbatch --array= ./sequencing/ATACSeq/preprocessing/calcENCODEQCMetrics.sh <sampleName>                  ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script calculates ENCODE library complexity metrics according to                                 ||
## https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#                              ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${METADIR}/samples.txt that lists sample names.                                                          ||
## - Variables in config file: ALIGNEDDIR                                                                             ||
## - Software: samtools, samstats, bedtools                                                                           ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> Name of sample specified in command line. Reads should be previously aligned                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           || 
## ${sampleName}.flagstat.qc, ${sampleName}.pbc.qc                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

sampleName=$1
echo "Calculating ENCODE QC metrics for" $sampleName
echo Job started on:
date -u

cd ${ALIGNEDDIR}

## If not found, create directory for ENCODEMetrics results
mkdir -p ENCODEMetrics

# run samstats on filtered, dedupped bam file if results file not found
if [ ! -s ENCODEMetrics/${sampleName}.flagstat.qc ]
then
  samtools sort -n --threads 10 ${sampleName}.filt.nodup.bam -O SAM  | SAMstats --sorted_sam_file -  --outf ENCODEMetrics/${sampleName}.flagstat.qc
  
else
  { echo "BAM file sorted, so not sorting"; exit 1;}
fi 

##  ================  ##
##    ENCODEMetrics   ##
##  ================  ##

if [ ! -s ENCODEMetrics/${sampleName}.pbc.qc ]
then
  
  # Results shown in file will be: TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
  
  #sort aligned reads by name
  samtools sort -n ${sampleName}.filt.dupmark.bam -o ${sampleName}.srt.tmp.bam
  
  # convert to bedPE and obtain fragment coordinates
  # sort by position and strand
  # obtain unique count statistics
  bedtools bamtobed -bedpe -i ${sampleName}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c \
  	| awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{if (mt > 0) {printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2} else printf "%d\t%d\t%d\t%d\t%s\t%s\t%s\n",mt,m0,m1,m2,NA,NA,NA }'  > ENCODEMetrics/${sampleName}.pbc.qc
  
  rm ${sampleName}.srt.tmp.bam
else
  { echo "Library complexity metrics calculated so not calculating"; exit 1;}
fi
