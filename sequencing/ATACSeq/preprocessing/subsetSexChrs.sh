## Performs peak calling for sex chr using shifted tagAlign files 

## EXECUTION
# sh ./ATACSeq/preprocessing/subsetSexChrs.sh <sampleName>
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR

## REQUIRES the following software
# samtools, gatk

## INPUT
# shifted tag align file
#

## OUTPUT
# ${ALIGNEDDIR}/sexChr/chrX.tn5.tagAlign.gz
# ${ALIGNEDDIR}/sexChr/chrY.tn5.tagAlign.gz
# ${ALIGNEDDIR}/SNPs/${sampleName}_chrX.gvcf


sampleName=$1

cd ${ALIGNEDDIR}/

mkdir -p SNPs

## variant calling by sample on X chromosome
## filter bam file to x chromosome
samtools view -b baseRecalibrate/${sampleName}_baserecal.bam chrX > baseRecalibrate/${sampleName}_baserecal_chrX.bam

samtools index baseRecalibrate/${sampleName}_baserecal_chrX.bam

gatk HaplotypeCaller \
-I baseRecalibrate/${sampleName}_baserecal_chrX.bam \
-R ${GENOMEFASTA} \
-O SNPs/${sampleName}_chrX.gvcf \
-ERC GVCF 

rm baseRecalibrate/${sampleName}_baserecal_chrX*

# for peak calling on sex chr, subset tagalign files
tagalign=${sampleName}.tn5.tagAlign.gz

mkdir -p sexChr

## extract X chr from tag align file
zcat $tagalign | grep "chrX"  | gzip -c > sexChr/${tagalign/tn5/chrX.tn5}

## extract Y chr from tag align file
zcat $tagalign | grep "chrY"  | gzip -c > sexChr/${tagalign/tn5/chrY.tn5}



