## collate flagstat output

## EXECUTION
# sh ChIPSeq/collateFlagStatOutput.sh 
# assumes config file has been loaded
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR

## REQUIRES the following software
# -

## INPUT
# -

## OUTPUT
# ${ALIGNEDDIR}/ENCODEMetrics/collateFlagStatMetrics.txt


cd ${ALIGNEDDIR}/ENCODEMetrics/

for file in *flagstat.qc;
do
	echo -n ${file%.flagstat.qc},
   awk 'BEGIN { ORS = "," } {print $1} END { printf( "\n" ); }' ${file}
   
done  > collateFlagStatMetrics.txt

