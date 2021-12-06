## for samples with multiple sequencing runs merge fastq files prior to processing
## move original fastq files into separate folder so that ignored through processing

mkdir -p ${DATADIRPE}/FastqParts

while read line
do 
	f1Merge=($(find ${DATADIRPE}/*/01_raw_reads/ -name *${line}*[rR]1*q.gz))
	## some files have a _b in filename
	f1Merge+=($(find ${DATADIRPE}/*/01_raw_reads/ -name *${line/_1_/_1_b_}*[rR]1*q.gz))
	FOLDER=$(dirname ${f1Merge[0]})
	echo ${f1Merge[@]}
	cat ${f1Merge[@]} > ${FOLDER}/${line}_r1.fq.gz
	
	mv ${f1Merge[@]} ${DATADIRPE}/FastqParts

	f2Merge=($(find ${DATADIRPE}/*/01_raw_reads/ -name *${line}*[rR]2*q.gz))
	f2Merge+=($(find ${DATADIRPE}/*/01_raw_reads/ -name *${line/_1_/_1_b_}*[rR]2*q.gz))
	FOLDER=$(dirname ${f2Merge[0]})
	echo ${f2Merge[@]}
	cat ${f2Merge[@]} > ${FOLDER}/${line}_r2.fq.gz
	
	mv ${f2Merge[@]} ${DATADIRPE}/FastqParts
	
done < ${DATADIRPE}/SamplesWithMultipleFastq.txt
