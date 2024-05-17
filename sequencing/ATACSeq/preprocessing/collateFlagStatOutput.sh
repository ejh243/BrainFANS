## =============================================================================================================##
##                          ATAC-seq pipeline STEP 4.2: Collate flagstat output                                 ##
## =============================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/collateFlagStatOutput.sh                                    ||
## - execute from scripts directory                                                                             ||
##                                                                                                              ||
## DESCRIPTION: This script collates the results from flagstat output                                           ||
##                                                                                                              ||
## OUTPUTS:                                                                                                     ||
## - collateFlagStatMetrics.txt                                                                                 ||
##                                                                                                              ||
## REQUIRES:                                                                                                    ||
## - bedtools, samtools                                                                                         ||
## - ALIGNEDDIR, PEAKDIR                                                                                        ||
##                                                                                                              ||
## =============================================================================================================##

cd ${ALIGNEDDIR}/ENCODEMetrics/

for file in *flagstat.qc;
do
	echo -n ${file%.flagstat.qc},
   awk 'BEGIN { ORS = "," } {print $1} END { printf( "\n" ); }' ${file}
   
done  > collateFlagStatMetrics.txt

