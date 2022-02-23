## Counts reads in GENCODE genes


## EXECUTION
# sh ./ATACSeq/preprocessing/geneBodyCounts.sh
# where 

# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR,PEAKCOUNTDIR,

## REQUIRES the following software
# featureCounts 

## INPUT
# indexed, bam files 

## OUTPUT
# ${PEAKCOUNTDIR}/GENCODE/geneCounts
# ${PEAKCOUNTDIR}/GENCODE/geneCounts.summary

cd ${ALIGNEDDIR}/


PULLDOWN=$(ls *PC*.GRCh38.karyo.deduplicated.bam)
CONTROL=$(ls *IC*.GRCh38.karyo.deduplicated.bam)

## count reads in genes

mkdir -p ${PEAKCOUNTDIR}/GENCODE/
featureCounts -T 8 -p -B -O -a ${GENCODEGTF} -o ${PEAKCOUNTDIR}/GENCODE/geneCounts ${PULLDOWN} ${CONTROL}




