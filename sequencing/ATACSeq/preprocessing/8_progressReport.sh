## searches for output files from each part of the QC pipeline and reports number of files found
## then searches specifically for each sample to identify what is missing

## EXECUTION
# sh ./ATACSeq/preprocessing/8_progressReport.sh
# assumes config file has been loaded
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# RAWDATADIR, FASTQCDIR, ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# -

## INPUT
# -

## OUTPUT
# ${METADIR}/SummariseSampleProcessingProgress.csv


echo "Running progress report"

## create array of all fastq files
cd ${RAWDATADIR}
FQFILES=($(find . -name '*[rR]1*q.gz')) ## this command searches for all fq files within

### First count number of output files for each part of the pipeline

echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""    

## check for FASTQC output
## NB two files per sample
echo "Number of fastqc reports found " $(ls ${FASTQCDIR}/*_fastqc.zip | wc -l)

## check for FASTP output
echo "Number of fastp reports found " $(ls ${TRIMDIR}/fastp_reports/*.json | wc -l)

## check for bowtie output
echo "Number of Bowtie reports found " $(ls ${ALIGNEDDIR}/*.bowtie.log | wc -l)

echo "Number of filtered aligned files found " $(ls ${ALIGNEDDIR}/*_postFilter_statsperchr.txt | wc -l)

## check for output of ENCODE QC metrics calculation
## first delete empty files 

find ${ALIGNEDDIR}/ENCODEMetrics/ -size  0 -print -delete

echo "Number of ENCODE QC metric output files found " $(ls ${ALIGNEDDIR}/ENCODEMetrics/*chr1*.pbc.qc | wc -l)


## check for peak calling output

echo "Number of MACS2 peak files (shifted tag align) found " $(ls ${PEAKDIR}/MACS/ShiftedTagAlign/*_peaks.broadPeak | wc -l)

echo "Number of MACS2 peak files (paired end) found " $(ls ${PEAKDIR}/MACS/BAMPE/*_peaks.broadPeak | wc -l)


## check for individual samples
## save output in txt file
echo "SampleID,DataFolder,R1Filename,R2Filename,FASTQCR1,FASTQCR2,FASTP,BOWTIE,FilteredAligned,ENCODEMetrics,MACS2Peaks1,MACS2Peaks2" > ${METADIR}/SummariseSampleProcessingProgress.csv
for f1 in ${FQFILES[@]}
do 
    sampleName=$(basename ${f1%_[rR][12]*}) ## sample name is everything before either r1 or R1
    ## later samples have an additional _S[num] in the file name need to remove
    sampleName=${sampleName%_S[0-9]*}
    echo "Processing" ${sampleName}
    
	## create filename for paired fastq file
    f2=$(basename $(ls ${RAWDATADIR}/${sampleName}*[rR]2*q.gz))    
    echo -n ${sampleName},${RAWDATADIR},$(basename ${f1}),${f2}, >> ${METADIR}/SummariseSampleProcessingProgress.csv
    
    if [ ! -s ${FASTQCDIR}/${sampleName}*[rR]1*fastqc.zip ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${FASTQCDIR}/${sampleName}*[rR]2*fastqc.zip ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi

    if [ ! -s ${TRIMDIR}/fastp_reports/${sampleName}*.json ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNEDDIR}/${sampleName}.bowtie.log ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${ALIGNEDDIR}/ENCODEMetrics/${sampleName}*chr1*.pbc.qc ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak.filt ]
    then
        echo -n "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo -n "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
    
    if [ ! -s ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak.filt ]
    then
        echo "N," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    else
        echo "Y," >> ${METADIR}/SummariseSampleProcessingProgress.csv
    fi
  
done
