## Written by Eilis
## Script to perform ${CNVDIR} calling with Penn${CNVDIR}
## generally following user guide: https://github.com/WGLab/Penn${CNVDIR}/blob/master/docs/user-guide/input.md

## assumes executed from script directory
mkdir -p ${CNVDIR}/PennCNVInput/

## create input files from final report
## replace array ids with sample IDs nb change "/" within names to "_"
${PENNCNVDIR}/split_illumina_report.pl -prefix ${CNVDIR}/PennCNVInput/ -revised_file ${CNVDIR}/ID_Map.txt ${RAWDATADIR}/${FILEPREFIX}FinalReport.txt

## need to replace column headers
FILES=${CNVDIR}/PennCNVInput/*
for f in $FILES
do
	sed -i '1s/.*/Name\tLog R Ratio\tB Allele Freq/' ${f}
done

##create snp position file excluding markers without a position
awk -v OFS='\t' '{if ($1 != 0) print $2, $1,$4}' ${RAWDATADIR}/${FILEPREFIX}.bim > ${PENNCNVDIR}/lib/gsa.gr38.snpposfile.tmp
echo -e "Name\tChr\tPos" | cat - ${PENNCNVDIR}/lib/gsa.gr38.snpposfile.tmp > ${PENNCNVDIR}/lib/gsa.gr38.snpposfile
rm ${PENNCNVDIR}/lib/gsa.gr38.snpposfile.tmp


## create pfb files
cd ${PENNCNVDIR}
mkdir -p ./pfb
cd ${CNVDIR}/PennCNVInput/

## instead use list of file 
${PENNCNVDIR}/compile_pfb.pl -listfile ${CNVDIR}/Samples.txt --snpposfile ${PENNCNVDIR}/lib/gsa.gr38.snpposfile -output ${PENNCNVDIR}/pfb/gsa.gr38.pfb


## create gc model file
cd ${PENNCNVDIR}/gc_file/
#gunzip hg38.gc5Base.txt.gz
sort -k 2,2 -k 3,3n hg38.gc5Base.txt > hg38.gc5BaseSorted.txt
cd ${CNVDIR}/PennCNVInput
${PENNCNVDIR}/cal_gc_snp.pl ${PENNCNVDIR}/gc_file/hg38.gc5BaseSorted.txt ${PENNCNVDIR}/lib/gsa.gr38.snpposfile -output ${PENNCNVDIR}/lib/GSA.gcmodel;


## call ${CNVDIR}s
cd ${CNVDIR}/
mkdir -p PennCNVOutput

## only run ${CNVDIR}s on samples that pass SNP QC
cd PennCNVInput/
${PENNCNVDIR}/detect_cnv.pl --test -hmm ${PENNCNVDIR}/lib/hh550.hmm -pfb ${PENNCNVDIR}/pfb/gsa.gr38.pfb --listfile ${CNVDIR}/Samples.txt -log ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_Merged.log -out ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_GCModel.rawcnv -gcmodel ${PENNCNVDIR}/lib/GSA.gcmodel;

## generate sample QC report 
${PENNCNVDIR}/filter_cnv.pl --qcsumout ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_SampleQC.txt --qclogfile ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_Merged.log ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_GCModel.rawcnv




