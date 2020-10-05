## Written by Eilis
## Uses shifted tagAlign files
## calls peaks across all samples for each fraction
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


mkdir -p ${PEAKDIR}
FRACTIONS=(N S T)

for type in ${FRACTIONS[@]};
do
	TAGFILES=$(ls ${ALIGNEDDIR}/*.tn5.tagAlign.gz | grep "ATAC-"${type})
	macs2 callpeak -t ${TAGFILES[@]} --outdir ${PEAKDIR} -n ${type} -f BED -g 2.9e9  -B --keep-dup all --shift 100 --extsize 200 --nomodel --broad

done

## run in all samples for unified peak set
TAGFILES=$(ls ${ALIGNEDDIR}/*.tn5.tagAlign.gz)
macs2 callpeak -t ${TAGFILES[@]} --outdir ${PEAKDIR} -n AllFractions -f BED -g 2.9e9  -B --keep-dup all --shift 100 --extsize 200 --nomodel --broad


## run MACS2 differential analysis to compare peaks between cell fractions

d1=$(grep "total tags in treatment" ${PEAKDIR}/N_peaks.xls | cut --delim=" " -f 6)
d2=$(grep "total tags in treatment" ${PEAKDIR}/S_peaks.xls | cut --delim=" " -f 6)
d3=$(grep "total tags in treatment" ${PEAKDIR}/T_peaks.xls | cut --delim=" " -f 6)
d1=`expr $d1 / 1000000`
d2=`expr $d2 / 1000000`
d3=`expr $d3 / 1000000`


macs2 bdgdiff --t1 ${PEAKDIR}/N_treat_pileup.bdg --c1 ${PEAKDIR}/N_control_lambda.bdg --t2 ${PEAKDIR}/S_treat_pileup.bdg\
   --c2 ${PEAKDIR}/S_control_lambda.bdg --d1 ${d1} --d2 ${d2} -g 100 -l 200 --o-prefix neun_vs_sox
   
macs2 bdgdiff --t1 ${PEAKDIR}/N_treat_pileup.bdg --c1 ${PEAKDIR}/N_control_lambda.bdg --t2 ${PEAKDIR}/T_treat_pileup.bdg\
   --c2 ${PEAKDIR}/T_control_lambda.bdg --d1 ${d1} --d2 ${d3} -g 100 -l 200 --o-prefix neun_vs_tot
   
macs2 bdgdiff --t1 ${PEAKDIR}/T_treat_pileup.bdg --c1 ${PEAKDIR}/T_control_lambda.bdg --t2 ${PEAKDIR}/S_treat_pileup.bdg\
   --c2 ${PEAKDIR}/S_control_lambda.bdg --d1 ${d3} --d2 ${d2} -g 100 -l 200 --o-prefix tot_vs_sox

