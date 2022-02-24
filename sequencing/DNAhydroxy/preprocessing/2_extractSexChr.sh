## Extract reads from X and Y chromosomes to perform sex checks


## EXECUTION
# sh ./ATACSeq/preprocessing/2_extractSexChr.sh <sampleName>
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, 

## REQUIRES the following software
# samtools

## INPUT
# aligned, filtered bam files for enriched and control fraction *.GRCh38.karyo.deduplicated.bam


## OUTPUT
# indexed bam files with just X and Y chromosomes

sampleName=$1

cd ${ALIGNEDDIR}

mkdir -p ${ALIGNEDDIR}/sexChr 

PCfile=$(ls ${sampleName}*PC*.GRCh38.karyo.deduplicated.bam)
ICfile=$(ls ${sampleName}*IC*.GRCh38.karyo.deduplicated.bam)

## read counts on sex chromosomes
samtools idxstats ${PCfile} | grep -w "chr[X|Y|M]" > ${ALIGNEDDIR}/sexChr/chrReadCounts_${sampleName}_PulldownBam.txt
samtools idxstats ${ICfile} | grep -w "chr[X|Y|M]"  > ${ALIGNEDDIR}/sexChr/chrReadCounts_${sampleName}_ControlBam.txt



## extract sex chromosomes for peak calling
samtools view -b ${PCfile} chrX chrY > ${ALIGNEDDIR}/sexChr/${PCfile}
samtools view -b ${ICfile} chrX chrY > ${ALIGNEDDIR}/sexChr/${ICfile}

## index
samtools index ${ALIGNEDDIR}/sexChr/${PCfile}
samtools index ${ALIGNEDDIR}/sexChr/${ICfile}