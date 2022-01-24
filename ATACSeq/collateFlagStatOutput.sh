## collate flagstat output

cd ${ALIGNEDDIR}/ENCODEMetrics/

for file in *flagstat.qc;
do
	echo -n ${file%_sorted_*flagstat.qc},
   awk 'BEGIN { ORS = "," } {print $1} END { printf( "\n" ); }' ${file}
   
done  > CollateFlagStatMetrics.txt

