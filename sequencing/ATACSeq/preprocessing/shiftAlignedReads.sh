## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 3.1: Shift Reads                                                ##
## ===================================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/shiftAlignedReads.sh <sampleName>                                 ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script converts bam files to tagalign files, calculates CC scores and shifts reads to perform    || 
## peak calling using single-end mode.  Adapted from: https://www.encodeproject.org/pipelines/ENCPL792NWO/            ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> sample name to shift reads and calculate CC scores                                              ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## ALIGNEDIR/*.tn5.tagAlign.gz, *.cc.qc, *.cc.plot.pdf                                                                ||
## - output will be peaks called in autosomal chromosomes and not aligned to blacklist regions                        ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for                                                                     || 
## - R run_spp.R script, bedtools, samtools                                                                           ||
## - ALIGNEDDIR                                                                                                       ||
##                                                                                                                    ||
## ===================================================================================================================##


sampleName=$1 
NTHREADS=8 
NREADS=15000000

cd ${ALIGNEDDIR}


# ===================
# Create tagAlign file
# ===================
#
# Create virtual SE file containing both read pairs
bedtools bamtobed -i ${sampleName}.filt.nodup.bam | awk 'BEGIN{OFS="\t"}{$5="1000";print $0}' | gzip -c > ${sampleName}.PE2SE.tn5.tagAlign.gz

# ================ 
# Create BEDPE file 
# ================ 
## requires name sorted bam file
samtools sort -n ${sampleName}.filt.nodup.bam > ${sampleName}.filt.nodup.nmsrt.bam
bedtools bamtobed -bedpe -mate1 -i ${sampleName}.filt.nodup.nmsrt.bam | gzip -c > ${sampleName}.filt.nodup.nmsrt.bedpe.gz
rm ${sampleName}.filt.nodup.nmsrt.bam

# =================================
# Subsample tagAlign file
# Restrict to one read end per pair for CC analysis
# ================================
zcat ${sampleName}.filt.nodup.nmsrt.bedpe.gz | grep -v “chrM” | shuf -n ${NREADS} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,"1000",$9}' | \
	gzip -c > "${sampleName}.filt.nodup.sample.$((NREADS /1000000)).MATE1.tagAlign.gz"

# =================================
# Calculate cross correlation scores
# ================================
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
Rscript ${PHANTOMPEAK}/run_spp.R -c="${sampleName}.filt.nodup.sample.$((NREADS /1000000)).MATE1.tagAlign.gz" -p=${NTHREADS} \
	-filtchr=chrM -savp=${sampleName}.subsample.cc.plot.pdf -out=${sampleName}.subsample.cc.qc -rf
sed -r 's/,[^\t]+//g' ${sampleName}.subsample.cc.qc > temp
mv temp ${sampleName}.subsample.cc.qc

echo 'Calculated cross correlation'

rm "${sampleName}.filt.nodup.sample.$((NREADS /1000000)).MATE1.tagAlign.gz"
rm ${sampleName}.filt.nodup.nmsrt.bedpe.gz

# ================
# Shift tagAlign file
# ================

zcat ${sampleName}.PE2SE.tn5.tagAlign.gz | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${sampleName}.tn5.tagAlign.gz

rm ${sampleName}.PE2SE.tn5.tagAlign.gz


