## Written by Eilis
## QC CNV calls including filtering and merging


# create plots of QC metrics and calculate exclusion thresholds (nb outliers defined as > 3 SD from mean)
QCFILTER=($(Rscript utilitys/summarizeCNVMetrics.r ${CNVDIR}/PennCNVOutput/SampleQC.txt ${PROCESSDIR}/PredictedPopulations.csv))

cd ${CNVDIR}/


# Filter samples
# Samples are excluded if > 10mbp of CNVs and outliers for Log2 ratio standard deviation, B-allele frequency drift, wave factor and total number of CNVs.
${PENNCNVPATH}/filter_cnv.pl --qclogfile PennCNVOutput/${FILEPREFIX}_Merged.log --qclrrsd ${QCFILTER[0]} --qcbafdrift ${QCFILTER[1]}  --qcwf ${QCFILTER[2]} --qcnumcnv ${QCFILTER[3]} --maxtotalcnvlength 10m --qcpassout PennCNVOutput/SamplesIncluded.txt --output PennCNVOutput/${FILEPREFIX}_GCModel_SampleFiltered.rawcnv PennCNVOutput/${FILEPREFIX}_GCModel.rawcnv 

# merge adjacent CNVs if separated by <25% of their combined length
${PENNCNVPATH}/clean_cnv.pl combineseg PennCNVOutput/${FILEPREFIX}_GCModel_SampleFiltered.rawcnv --signalfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile --output PennCNVOutput/${FILEPREFIX}_GCModel_Merged.rawcnv --fraction 0.25

# Filter CNVs
# CNVs are excluded if covered by less than 10 probes, less than 15kb in length
${PENNCNVPATH}/filter_cnv.pl --numsnp 10 --length 15k --output PennCNVOutput/${FILEPREFIX}_GCModel_MergedFiltered.rawcnv PennCNVOutput/${FILEPREFIX}_GCModel_Merged.rawcnv 


${PENNCNVPATH}/scan_region.pl PennCNVOutput/${FILEPREFIX}_GCModel_MergedFiltered.rawcnv --knowngene ${UCSCUTILS}/knownGene.Gencodev39.hg38 --kgxref ${UCSCUTILS}/kgXref.hg38 > PennCNVOutput/${FILEPREFIX}_GCModel_MergedFiltered_AnnoGencodev39.rawcnv

## Note needs to be done separately for each sample
#./visualize_cnv.pl -format plot -snpposfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile PennCNVOutput/${FILEPREFIX}_GCModel_MergedFiltered_AnnoGencodev29.rawcnv --signal 

