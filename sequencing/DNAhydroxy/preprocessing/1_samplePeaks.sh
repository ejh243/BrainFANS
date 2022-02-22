## Performs sample level peak calling using MACS and EPIC2/
## filters peaks to exclude black listed regions


## EXECUTION
# sh ./ATACSeq/preprocessing/1_samplePeaks.sh <sampleName>
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2, epic2

## INPUT
# aligned, filtered bam files for enriched and control fraction *.GRCh38.karyo.deduplicated.bam


## OUTPUT


sampleName=$1

cd ${ALIGNEDDIR}

mkdir -p ${PEAKDIR}/samplePeaks/MACS/

PCfile=$(ls ${sampleName}*PC*.GRCh38.karyo.deduplicated.bam)
ICfile=$(ls ${sampleName}*IC*.GRCh38.karyo.deduplicated.bam)	
macs2 callpeak -t ${PCfile} -c ${ICfile} -f BAM -g hs -n ${sampleName} -q 0.01 --outdir ${PEAKDIR}/samplePeaks/MACS/


mkdir -p ${PEAKDIR}/samplePeaks/EPIC2/

epic2 -t ${PCfile} -c ${ICfile} --genome hg38 -fdr 0.01 -o ${PEAKDIR}/samplePeaks/EPIC2/${sampleName}.txt


## exclude peaks aligned to blacklist regions
mkdir -p ${PEAKDIR}/samplePeaks/MACS/BlacklistFiltered/

bedtools intersect -a ${PEAKDIR}/samplePeaks/MACS/${sampleName}_peaks.narrowPeak -b ${BLACKLIST} -v > ${PEAKDIR}/samplePeaks/MACS/BlacklistFiltered/${sampleName}_peaks.narrowPeak.filt

rm ${PEAKDIR}/samplePeaks/MACS/${sampleName}*
 
mkdir -p ${PEAKDIR}/samplePeaks/EPIC2/BlacklistFiltered/
bedtools intersect -a ${PEAKDIR}/samplePeaks/EPIC2/${sampleName}.txt -b ${BLACKLIST} -v > ${PEAKDIR}/samplePeaks/EPIC2/BlacklistFiltered/${sampleName}_peaks.broadPeak.filt
  
rm ${PEAKDIR}/samplePeaks/EPIC2/${sampleName}*

## calc FRIP
mkdir -p  ${PEAKDIR}/QCoutput/
## save output in txt file
echo "SampleName,BAMPulldownReads,BAMControlReads,MACS2Peaks,PulldownReadsinMACS2Peaks,ControlReadsinMACS2Peaks,EPIC2Peaks,PulldownReadsinEPIC2Peaks,ControlReadsinEPIC2Peaks" > ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv

echo -n ${sampleName},  >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
# count total reads in peak input files
echo -n $(samtools view -c ${ALIGNEDDIR}/${PCfile}), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
echo -n $(samtools view -c ${ALIGNEDDIR}/${ICfile}), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv

# count reads in peak regions
## MACS2 peaks 
if [ -s ${PEAKDIR}/samplePeaks/MACS/BlacklistFiltered/${sampleName}_peaks.narrowPeak.filt ]
then
	echo -n $(wc -l ${PEAKDIR}/samplePeaks/MACS/BlacklistFiltered/${sampleName}_peaks.narrowPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
	echo -n $(bedtools sort -i ${PEAKDIR}/samplePeaks/MACS/BlacklistFiltered/${sampleName}_peaks.narrowPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${PCfile} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
	echo -n $(bedtools sort -i ${PEAKDIR}/samplePeaks/MACS/BlacklistFiltered/${sampleName}_peaks.narrowPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ICfile} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
else
	echo -n NA,NA,NA, >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
fi


## EPIC2 peaks 
if [ -s ${PEAKDIR}/samplePeaks/EPIC2/BlacklistFiltered/${sampleName}_peaks.broadPeak.filt ]
then
	echo -n $(wc -l ${PEAKDIR}/samplePeaks/EPIC2/BlacklistFiltered/${sampleName}_peaks.broadPeak.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
	echo -n $(bedtools sort -i ${PEAKDIR}/samplePeaks/EPIC2/BlacklistFiltered/${sampleName}_peaks.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${PCfile} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
	echo $(bedtools sort -i ${PEAKDIR}/samplePeaks/EPIC2/BlacklistFiltered/${sampleName}_peaks.broadPeak.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ICfile} -b stdin | wc -l | cut -f1 -d' ') >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
else
	echo NA,NA,NA >> ${PEAKDIR}/QCoutput/FRIP_${sampleName}.csv
fi
