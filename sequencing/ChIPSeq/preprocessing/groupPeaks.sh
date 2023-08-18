## Performs fraction/group level peak calling using MACS 
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./ChIPSeq/preprocessing/fractionPeaks.sh 

# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2

## INPUT
# array of bam files

## OUTPUT
# ${PEAKDIR}/MACS/*.broadPeak.filt

## broad and narrow marks taken from https://www.encodeproject.org/chip-seq/histone/
broad=(H3F3A H3K27me3 H3K36me3 H3K4me1 H3K79me2 H3K79me3 H3K9me1 H3K9me2 H4K20me1)
narrow=(H2AFZ H3ac H3K27ac H3K4me2 H3K4me3 H3K9ac)


cd ${ALIGNEDDIR}

SAMPLEFILE=$1
mark=$2

echo
echo 'Starting peak calling at:'
date -u

QC1SAMPLES=($(grep "${GROUP}" ${SAMPLEFILE} | awk '{print $1}' ))
CONTROL=($(grep "${GROUP}" ${SAMPLEFILE} | awk '{print $4}'))

FILES="${QC1SAMPLES[@]/%/*nodup.bam}"
CONTROL="${CONTROL[@]/%/*nodup.bam}"

echo 'Files for peak-calling are:' ${FILES}

if [[ ${broad[*]} =~  ${mark} ]]
then
	# call for broad peaks
	macs2 callpeak -t ${FILES} -c ${CONTROL} --outdir ${PEAKGROUP} -n ${GROUP}_${mark} -g 2.9e9 -q 1e-2 --broad

	## exclude peaks aligned to blacklist regions
	bedtools intersect -v -a ${PEAKGROUP}/${GROUP}_${mark}_peaks.broadPeak -b ${BLACKLIST} \
		| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
		| grep -P 'chr[\dXY]+[ \t]' > ${PEAKGROUP}/${GROUP}_${mark}.broadPeak.filt
elif [[ ${narrow[*]} =~  ${mark} ]]
then
  	# calculate narrowPeak
  	macs2 callpeak -t ${FILES} -c ${CONTROL} --outdir ${PEAKGROUP} -n ${GROUP}_${mark} -g 2.9e9 -q 1e-2 

  	## exclude peaks aligned to blacklist regions
  	bedtools intersect -v -a ${PEAKGROUP}/${GROUP}_${mark}_peaks.narrowPeak -b ${BLACKLIST} \
    	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    	| grep -P 'chr[\dXY]+[ \t]' > ${PEAKGROUP}/${GROUP}_${mark}.narrowPeak.filt

else
  	echo 'Target not found in either broad or narrow peak mark'
fi


if [[ $? == 0 ]]
then 
	echo 'Finished calling peaks in common'
fi
