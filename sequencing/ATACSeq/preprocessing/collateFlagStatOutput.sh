## collate flagstat output

## EXECUTION
# sh ATACSeq/collateFlagStatOutput.sh 
# assumes config file has been loaded
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# bedtools, samtools

## INPUT
# blacklist filtered peak lists
#

## OUTPUT
# ${PEAKDIR}/QCOutput/FRIP_*.csv


cd ${ALIGNEDDIR}/ENCODEMetrics/

for file in *flagstat.qc;
do
	echo -n ${file%.flagstat.qc},
   awk 'BEGIN { ORS = "," } {print $1} END { printf( "\n" ); }' ${file}
   
done  > collateFlagStatMetrics.txt

