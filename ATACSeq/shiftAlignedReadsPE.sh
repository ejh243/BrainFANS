

## this script takes filtered bam file converts to a tagalign file, calculate CC scores and shifts reads ready for peak calling
## based on code developed as part of the ENCODE ATAC-seq pipeline: https://www.encodeproject.org/pipelines/ENCPL792NWO/

# input Filtered BAM file and name sorted BAM file
NTHREADS=8 


cd ${ALIGNEDDIR}

OFPREFIX=$1
FINAL_BAM_FILE=${OFPREFIX}_depDup_q30.bam
FINAL_NMSRT_BAM_PREFIX=${OFPREFIX}.filt.nmsrt.nodup 
FINAL_NMSRT_BAM_FILE=${FINAL_NMSRT_BAM_PREFIX}.bam # To be stored 

# ===================
# Create tagAlign file
# ===================

# Create virtual SE file containing both read pairs
FINAL_TA_FILE="${OFPREFIX}.PE2SE.tagAlign.gz"
bedtools bamtobed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > ${FINAL_TA_FILE}

# ================ 
# Create BEDPE file 
# ================ 
FINAL_BEDPE_FILE="${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz" 

## requires name sorted bam file
samtools sort -n ${FINAL_BAM_FILE} > ${FINAL_NMSRT_BAM_FILE} 
bedtools bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | gzip -c > ${FINAL_BEDPE_FILE} 

# =================================
# Subsample tagAlign file
# Restrict to one read end per pair for CC analysis
# ================================
NREADS=15000000
SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.$((NREADS /1000000)).MATE1.tagAlign.gz"
zcat ${FINAL_BEDPE_FILE} | grep -v “chrM” | shuf -n ${NREADS} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"N","1000",$9}' | gzip -c > ${SUBSAMPLED_TA_FILE}

# =================================
# Calculate cross correlation scores
# ================================

CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"
# CC_SCORE FILE format
# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag
Rscript ${PHANTOMPEAK}run_spp.R -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE} -rf
sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > temp
mv temp ${CC_SCORES_FILE}


# ================
# Shift tagAlign file
# ================

shifted_tag=${OFPREFIX}.tn5.tagAlign.gz
zcat $FINAL_TA_FILE | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > ${shifted_tag}


