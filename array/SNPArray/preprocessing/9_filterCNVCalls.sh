## Written by Eilis
## QC CNV calls including filtering and merging


cd ${DATADIR}/scripts
CNVPATH=${DATADIR}/SNPdata/CNV


# create plots of QC metrics and calculate exclusion thresholds (nb outliers defined as > 3 SD from mean)
QCFILTER=($(Rscript SNPdata/summarizeCNVMetrics.r ${CNVPATH}/PennCNVOutput/SampleQC.txt ../SNPdata/Merged/PredictedPopulations.csv))

cd ${DATADIR}/SNPdata/CNV


# Filter samples
# Samples are excluded if > 10mbp of CNVs and outliers for Log2 ratio standard deviation, B-allele frequency drift, wave factor and total number of CNVs.
${PENNCNVPATH}/filter_cnv.pl --qclogfile PennCNVOutput/SCZ2.log --qclrrsd ${QCFILTER[0]} --qcbafdrift ${QCFILTER[1]}  --qcwf ${QCFILTER[2]} --qcnumcnv ${QCFILTER[3]} --maxtotalcnvlength 10m --qcpassout PennCNVOutput/SamplesIncluded.txt --output PennCNVOutput/SCZ_GCModel_SampleFiltered.rawcnv PennCNVOutput/SCZ_GCModel.rawcnv 

# merge adjacent CNVs if separated by <25% of their combined length
${PENNCNVPATH}/clean_cnv.pl combineseg PennCNVOutput/SCZ_GCModel_SampleFiltered.rawcnv --signalfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile --output PennCNVOutput/SCZ_GCModel_Merged.rawcnv --fraction 0.25

# Filter CNVs
# CNVs are excluded if covered by less than 10 probes, less than 15kb in length
${PENNCNVPATH}/filter_cnv.pl --numsnp 10 --length 15k --output PennCNVOutput/SCZ_GCModel_MergedFiltered.rawcnv PennCNVOutput/SCZ_GCModel_Merged.rawcnv 


${PENNCNVPATH}/scan_region.pl PennCNVOutput/SCZ_GCModel_MergedFiltered.rawcnv --knowngene ${DATADIR}/References/UCSC/knownGene.Gencodev29.hg38 --kgxref ${DATADIR}/References/UCSC/kgXref.hg38 > PennCNVOutput/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv

## Note needs to be done separately for each sample
#./visualize_cnv.pl -format plot -snpposfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile PennCNVOutput/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv --signal 


cd ${DATADIR}/scripts

QCFILTER=($(Rscript SNPdata/summarizeCNVMetrics.r ${CNVPATH}/PennCNVOutput/EURonly/SampleQC.txt ../SNPdata/Merged/PredictedPopulations.csv))


cd ${DATADIR}/SNPdata/CNV
# Filter samples
# Samples are excluded if > 10mbp of CNVs and outliers for Log2 ratio standard deviation, B-allele frequency drift, wave factor and total number of CNVs.
${PENNCNVPATH}/filter_cnv.pl --qclogfile PennCNVOutput/EURonly/SCZ2.log --qclrrsd ${QCFILTER[0]} --qcbafdrift ${QCFILTER[1]}  --qcwf ${QCFILTER[2]} --qcnumcnv ${QCFILTER[3]} --maxtotalcnvlength 10m --qcpassout PennCNVOutput/EURonly/SamplesIncluded.txt --output PennCNVOutput/EURonly/SCZ_GCModel_SampleFiltered.rawcnv PennCNVOutput/EURonly/SCZ_GCModel.rawcnv 

# merge adjacent CNVs if separated by <25% of their combined length
${PENNCNVPATH}/clean_cnv.pl combineseg PennCNVOutput/EURonly/SCZ_GCModel_SampleFiltered.rawcnv --signalfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile --output PennCNVOutput/EURonly/SCZ_GCModel_Merged.rawcnv --fraction 0.25

# Filter CNVs
# CNVs are excluded if covered by less than 10 probes, less than 15kb in length
${PENNCNVPATH}/filter_cnv.pl --numsnp 10 --length 15k --output PennCNVOutput/EURonly/SCZ_GCModel_MergedFiltered.rawcnv PennCNVOutput/EURonly/SCZ_GCModel_Merged.rawcnv 


${PENNCNVPATH}/scan_region.pl PennCNVOutput/EURonly/SCZ_GCModel_MergedFiltered.rawcnv --knowngene ${DATADIR}/References/UCSC/knownGene.Gencodev29.hg38 --kgxref ${DATADIR}/References/UCSC/kgXref.hg38 > PennCNVOutput/EURonly/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv
