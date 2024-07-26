#!/bin/bash

## =====================================================================================================##
##               ATAC-seq pipeline STEP 6.1: Verify Genotype                                            ##
## =====================================================================================================##
## EXECUTION: sh ./subScripts/compareBamWithGenotype.sh <sampleName>  <vcf id>                          ||
## - execute from pipeline's main scripts directory                                                     ||
##                                                                                                      ||
## DESCRIPTION: This script compares genotype information of sample with existing genotype SNPs info to ||
## verify there is no contamination.                                                                    ||
##                                                                                                      ||
## INPUTS:                                                                                              || 
## $1 -> <sampleName> sample name verify genotype with their matched VCF file                           ||
## $2 -> <vcf id> ID within the vcf file for the input sample                                           ||
##                                                                                                      ||
## OUTPUTS:                                                                                             ||
## - ALIGNED_DIR/baseRecalibrate/<sampleName>_metrics.txt,/<sampleName>_baserecal.bam,                  ||
##   /<sampleName>_baserecal.bai                                                                        ||
## - ALIGNED_DIR/genotypeConcordance/<sampleName>.selfSM, /<sampleName>.log, /<sampleName>.depthSM      ||
## - Samples for which a matching VCF is not found will be printed to: $META_DIR/noVCFfound.txt         ||
##                                                                                                      ||
## REQUIRES:                                                                                            ||
## - sorted bam file, VCF file with SNP chip data                                                       ||
## - Softwares: gatk, samtools, picard, verifyBamID                                                     ||
## - Variables in config file: ALIGNED_DIR, GENOMEFASTA, KGREF,GENODIR,VERIFYBAMID                      ||
## - Reference genome: GENOMEFASTA: genome.fa                                                           ||
##                                                                                                      ||
## =====================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

cd ${ALIGNED_DIR}/

## If directory does not exist, create
mkdir -p ${ALIGNED_DIR}/baseRecalibrate

sampleName=$1
vcfid=$2

## ============ ##
##    VERIFY    ##
## ============ ##

if [[ "${vcfid}" != "" ]]
then 
    echo "processing" ${sampleName} "with" ${vcfid}
    bamfile=${sampleName}_sorted.bam

    re='^[0-9]{5}'
    if [[ $sampleName =~ $re ]] ;
    then
       projectID=${sampleName%%_*}
    else
       projectID="NA"
    fi

    ## mark duplicates only
    java -jar $EBROOTPICARD/picard.jar  MarkDuplicates INPUT=${ALIGNED_DIR}/${bamfile} OUTPUT=${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup.bam METRICS_FILE=${ALIGNED_DIR}/baseRecalibrate/${sampleName}_metrics.txt     

    ## add read group
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
           I=${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup.bam \
           O=${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup_rglabelled.bam \
           RGID=${projectID} \
           RGLB=lib1 \
           RGPL=ILLUMINA \
           RGPU=unit1 \
           RGSM=${sampleName}

    rm ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup.bam

    ##index
    samtools index ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup_rglabelled.bam

    # recalibrate bases in bam files
    gatk BaseRecalibrator \
        -R ${GENOMEFASTA}/genome.fa \
        -I ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup_rglabelled.bam \
        --known-sites ${KGREF}/1000G_omni2.5.hg38.vcf.gz \
        -O ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_recal_data.table

    gatk ApplyBQSR \
      -R ${GENOMEFASTA}/genome.fa \
       -I ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup_rglabelled.bam \
       --bqsr-recal-file ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_recal_data.table \
       -O ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_baserecal.bam
       
    rm ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_dedup_rglabelled.bam*
	  rm ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_recal_data.table

    mkdir -p ${ALIGNED_DIR}/genotypeConcordance
    
    ${VERIFYBAMID} --vcf ${SNP} --bam ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_baserecal.bam --out ${ALIGNED_DIR}/genotypeConcordance/${sampleName} --ignoreRG --smID ${vcfid} --self
else
    echo "No genotype data available"
    echo $sampleName >> ${META_DIR}/noVCFfound.txt
fi

if [[ ! -f ${ALIGNED_DIR}/genotypeConcordance/${sampleName}.selfSM ]]
then
  { echo "Genotype could not verified for ${sampleName}. Check whether a matching VCF IDs exists for this sample." ; exit 1 ; }
else
  echo "Genotype verification done for ${sampleName}."
  echo Job finished at: 
  date -u
fi




 