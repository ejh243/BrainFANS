## run at relaxed q-value threshold to enable classification of null regions in the genome

type_input=$1


type_control=${type_input/Input/Control}
type=$(basename $type_input)
type=${type%.txt}
type=${type##PeakCallingInputFiles}
	
input_samples=$(cat ${type_input})
control_samples=$(cat ${type_control})

cd ${ALIGNEDDIR}
	
macs2 callpeak -t ${input_samples[@]} -c ${control_samples[@]} -f BAM -g hs -n ${type} -B -q 0.5 --outdir ${PEAKDIR}/MACS2/	

## filter out blacklist peaks

bedtools intersect -a ${PEAKDIR}/MACS2/${type}_peaks.narrowPeak -b ${REFDIR}/ENCODE/Blacklist/hg38-blacklist.v2.bed -v > ${PEAKDIR}/MACS2/${type}_peaks.narrowPeak.filt



	