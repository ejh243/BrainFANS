## Written by Eilis
## Uses HMMRATAC 
## calls peaks per sample for QC purposes
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331

cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}/HMMRATAC

## create file paths from sample name
## requires sorted bam file with index from paired end data
sampleName=$1
BAM=${sampleName}_depDup_q30.bam
INDEX=${sampleName}_depDup_q30.bam.bai


java -jar ${HMMRATAC} -b ${BAM} -i ${INDEX} -g ${GENOMESIZE} -e ${BLACKLIST} -o ${PEAKDIR}/HMMRATAC/${sampleName}