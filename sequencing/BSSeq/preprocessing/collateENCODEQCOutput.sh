## collate encode qc output

## EXECUTION
# sh 
# assumes config file has been loaded
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR


## OUTPUT
# ${ALIGNEDDIR}/ENCODEMetrics/collateFlagStatMetrics.txt


cd ${ALIGNEDDIR}/ENCODEMetrics/

echo 'sampleID,coverage,convLambda,convCHG'> collateENCODEQCMetrics.txt

for file in *.qc;
do
	sampleName=${file%.qc}
	var=$( grep 'C methylated in CHG context:' ../bismarkReports/${sampleName}*.txt | awk '{ print $6 }' | rev | cut -c2- | rev )
	echo -n ${sampleName},
   	awk 'BEGIN { ORS = "," } {print $1} ' ${file}
   	echo "100 - $var" | bc 
   
done  >> collateENCODEQCMetrics.txt


# END { printf( "\n" ); }