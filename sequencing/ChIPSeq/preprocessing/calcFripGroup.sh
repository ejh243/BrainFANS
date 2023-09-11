FAIL=$@


echo 'passed:' $@
echo 'FAIL:' ${FAIL[@]}

for sampleName in ${FAIL[@]}
do
	# formulate peak input filenames
	aligned=${sampleName}*nodup.bam

	## save output in txt file
	echo "sampleName,bamTotalReads,samplePeaks,readsInSamplePeaks,subsetPeaks,readsInSubsetPeaks" > ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv

	echo -n ${sampleName}, >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv

	#sampletotalreads
	echo -n $(zcat ${ALIGNEDDIR}/${aligned} | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#samplepeaks
	echo -n $(wc -l ${PEAKDIR}/${sampleName}*.filt | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#reads in peaks
	echo -n $(bedtools sort -i ${PEAKDIR}/${sampleName}*.filt | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${aligned} -b stdin | wc -l | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#subsetpeaks
	echo -n $(wc -l ${PEAKDIR}/subset/*.narrowPeak | cut -f1 -d' '), >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
	#reads in subset peaks
	echo -n $(bedtools sort -i ${PEAKDIR}/subset/*.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a ${ALIGNEDDIR}/${aligned} -b stdin | wc -l | cut -f1 -d' ') >> ${PEAKDIR}/QCOutput/subset/FRIP_${sampleName}.csv
done