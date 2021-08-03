## run at relaxed q-value threshold to enable classification of null regions in the genome

cd ${ALIGNEDDIR}

mkdir -p ${PEAKDIR}

sampleSets=($(ls ../sampleLists/PeakCallingInputFiles*))

for type_input in ${sampleSets[@]}
do
	type_control=${type_input/Input/Control}
	type=$(basename $type_input)
	type=${type%.txt}
	type=${type##PeakCallingInputFiles}
	
	input_samples=$(cat ${type_input})
	control_samples=$(cat ${type_control})
	
	macs2 callpeak -t ${input_samples[@]} -c ${control_samples[@]} -f BAM -g hs -n ${type} -B -q 0.5 --outdir ${PEAKDIR}	
	
done

## run peaks across all smaples

input_samples=($(ls *PC*_L00.bml.GRCh38.karyo.deduplicated.bam))
control_samples=($(ls *IC*_L00.bml.GRCh38.karyo.deduplicated.bam))

macs2 callpeak -t ${input_samples[@]} -c ${control_samples[@]} -f BAM -g hs -n ALL -B -q 0.5 --outdir ${PEAKDIR}
