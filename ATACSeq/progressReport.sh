## searches for output files from each part of the QC pipeline and reports number of files found
## then searches specifically for each sample to identify what is missing


## load config file provided on command line when submitting job
echo "Loading config file: "
echo $1
source ./$1 

echo "Running progress report"


cd ${DATADIRPE}
## find all folders with fastq files in
RAWDATADIR=($(find . -name '01_raw_reads'))

## create array of all fastq files
FQFILES=()
for FOLDER in ${RAWDATADIR[@]}
do 	
	if [ ${FOLDER} != FastqParts ]
	then 
		FQFILES+=($(find ${FOLDER} -name '*[rR]1*q.gz')) ## this command searches for all fq files within
	fi
done

### First count number of outpu files for each part of the pipeline

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

## check for FASTQC output
## NB two files per sample
echo "Number of fastqc reports found " $(ls ${FASTQCDIR}/*_fastqc.zip | wc -l)

## check for FASTP output
echo "Number of fastp reports found " $(ls ${FOLDERTRIM}/fastp_reports/*.json | wc -l)

## check for bowtie output

echo "Number of Bowtie reports found " $(ls ${ALIGNEDDIR}/*.bowtie.log | wc -l)

echo "Number of filtered aligned files found " $(ls ${ALIGNEDDIR}/*_postFilter_statsperchr.txt | wc -l)

## check for output of ENCODE QC metrics calculation
## first delete empty files 

find ${ALIGNEDDIR}/ENCODEMetrics/ -size  0 -print -delete

echo "Number of ENCODE QC metric output files found " $(ls ${ALIGNEDDIR}/ENCODEMetrics/*chr1*.pbc.qc | wc -l)


## check for individual samples
## save output in txt file
echo "SampleID,DataFolder,R1Filename,R2Filename,FASTQCR1,FASTQCR,FASTP,BOWTIE,FilteredAligned,ENCODEMetrics" > ${DATADIRPE}/SummariseSampleProcessingProgress.csv
for f1 in ${FQFILES[@]}
do 
	FOLDER=$(dirname ${f1})
	f1=$(basename $f1)
	sampleName=${f1%_[rR]*}
	## later samples have an additional _S[num] in the file name need to remove
	sampleName=${sampleName%_S[0-9]*}
	## create filename for paired fastq file
	f2=$(basename $(ls ${FOLDER}/${sampleName}*[rR]2*q.gz))
	
	echo -n ${sampleName},${FOLDER}, ${f1},${f2}, >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	
	if [ ! -s ${FASTQCDIR}/${sampleName}*[rR]1*fastqc.zip ]
	then
		echo -n "N," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	else
		echo -n "Y," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	fi
	
	if [ ! -s ${FASTQCDIR}/${sampleName}*[rR]2*fastqc.zip ]
	then
		echo -n "N," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	else
		echo -n "Y," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	fi

	if [ ! -s ${FOLDERTRIM}/fastp_reports/${sampleName}*.json ]
	then
		echo -n "N," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	else
		echo -n "Y," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	fi
	
	if [ ! -s ${ALIGNEDDIR}/${sampleName}.bowtie.log ]
	then
		echo -n "N," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	else
		echo -n "Y," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	fi
	
	if [ ! -s ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt ]
	then
		echo -n "N," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	else
		echo -n "Y," >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	fi
	
	if [ ! -s ${ALIGNEDDIR}/ENCODEMetrics/${sampleName}*chr1*.pbc.qc ]
	then
		echo "N" >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	else
		echo "Y" >> ${DATADIRPE}/SummariseSampleProcessingProgress.csv
	fi
done
