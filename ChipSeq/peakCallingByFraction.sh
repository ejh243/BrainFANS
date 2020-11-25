## Written by Eilis

FRACTIONS=(NEU SOX TOT NEG)

for type in ${FRACTIONS[@]};
do
	SUBSET=$(ls ${ALIGNEDDIR}/*_depDup_q30.bam| grep ${type})

	macs2 callpeak -t ${SUBSET} --outdir ${PEAKDIR} -n ${type} -g 2.9e9  -B --nomodel --shift 0
done

d1=$(grep "total tags in treatment" ${PEAKDIR}/Neon_peaks.xls | cut --delim=" " -f 6)
d2=$(grep "total tags in treatment" ${PEAKDIR}/Sox_peaks.xls | cut --delim=" " -f 6)
d3=$(grep "total tags in treatment" ${PEAKDIR}/Total_peaks.xls | cut --delim=" " -f 6)
d1=`expr $d1 / 1000000`
d2=`expr $d2 / 1000000`
d3=`expr $d3 / 1000000`


macs2 bdgdiff --t1 ${PEAKDIR}/Neon_treat_pileup.bdg --c1 ${PEAKDIR}/Neon_control_lambda.bdg --t2 ${PEAKDIR}/Sox_treat_pileup.bdg\
   --c2 ${PEAKDIR}/Sox_control_lambda.bdg --d1 ${d1} --d2 ${d2} -g 100 -l 200 --o-prefix neun_vs_sox
   
macs2 bdgdiff --t1 ${PEAKDIR}/Neon_treat_pileup.bdg --c1 ${PEAKDIR}/Neon_control_lambda.bdg --t2 ${PEAKDIR}/Total_treat_pileup.bdg\
   --c2 ${PEAKDIR}/Total_control_lambda.bdg --d1 ${d1} --d2 ${d3} -g 100 -l 200 --o-prefix neun_vs_tot
   
macs2 bdgdiff --t1 ${PEAKDIR}/Total_treat_pileup.bdg --c1 ${PEAKDIR}/Total_control_lambda.bdg --t2 ${PEAKDIR}/Sox_treat_pileup.bdg\
   --c2 ${PEAKDIR}/Sox_control_lambda.bdg --d1 ${d3} --d2 ${d2} -g 100 -l 200 --o-prefix tot_vs_sox


