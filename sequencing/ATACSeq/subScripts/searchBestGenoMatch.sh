#!/bin/bash
## =============================================================================================##
##                    ATAC-seq pipeline STEP 6.3: Best Genotype                                 ##
## =============================================================================================##
## EXECUTION: sh ./subScripts/searchBestGenotype.sh <sampleName>  <vcf id>                      ||
## - execute from pipeline's main scripts directory                                             ||
##                                                                                              ||
## DESCRIPTION: This script takes potential contaminated genotype samples and aims to find the  ||
## best matching genotype                                                                       ||
##                                                                                              ||
## INPUTS:                                                                                      || 
## $1 -> <sampleName> sample name verify genotype with their matched VCF file                   ||
## $2 -> <vcf id> ID within the vcf file for the input sample                                   ||
##                                                                                              ||                                                                                                                    
## OUTPUTS:                                                                                     ||
## - ALIGNED_DIR/genotypeConcordance/<sampleName>.selfSM, /<sampleName>.log,                    ||                
##   /<sampleName>.depthSM                                                                      ||                                        
## REQUIRES:                                                                                    ||
## - sorted bam file, VCF file with SNP chip data                                               ||
## - Softwares: verifyBamID                                                                     ||
## - Variables in config file: ALIGNED_DIR, VERIFYBAMID                                         ||
##                                                                                              ||
## =============================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

sampleName=$1
cd ${ALIGNED_DIR}/genotypeConcordance/
sampleName=${sampleName%.selfSM}
vcfid=$2
vcfid="${vcfid%"${vcfid##*[![:space:]]}"}" 

cd ${ALIGNED_DIR}/
echo "processing" ${sampleName} "with" ${vcfid}

## ======================== ##
##    FIND BEST GENOTYPE    ##
## ======================== ##

${VERIFYBAMID} --vcf ${SNP} --bam ${ALIGNED_DIR}/baseRecalibrate/${sampleName}_baserecal.bam --out ${ALIGNED_DIR}/genotypeConcordance/${sampleName} --ignoreRG --smID ${vcfid} --best

if [[ ! -f ${ALIGNED_DIR}/genotypeConcordance/${sampleName}.bestSM ]]
then
  { echo "The best Genotype could not be found for ${sampleName}. Check whether a matching VCF IDs exists for this sample." ; exit 1 ;}
else
  echo "Best genotype found for ${sampleName}."
  echo Job finished at: 
  date -u
fi
