## Written by Eilis
## Script to perform CNV calling with PennCNV
## generally following user guide: https://github.com/WGLab/PennCNV/blob/master/docs/user-guide/input.md

## assumes executed from script directory
cd ${DATADIR}/
cd SNPdata/

## create input files from final report
## replace array ids with sample IDs nb change "/" within names to "_"
${PENNCNVPATH}/split_illumina_report.pl -prefix CNV/PennCNVInput/ -revised_file CNV/IDMap.txt SCZ2_gr38_FinalReport.txt
## need to replace column headers
FILES=CNV/PennCNVInput/*
for f in $FILES
do
	sed -i '1s/.*/Name\tLog R Ratio\tB Allele Freq/' ${f}
done

##create snp position file excluding markers without a position
awk -v OFS='\t' '{if ($1 != 0) print $2, $1,$4}' SCZ2_gr38_binary.bim > ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp
echo -e "Name\tChr\tPos" | cat - ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp > ${PENNCNVPATH}/lib/gsa.gr38.snpposfile
rm ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp


## create pfb files
cd ${PENNCNVPATH}
mkdir -p ./pfb
cd ${DATADIR}/SNPdata/CNV/PennCNVInput/

## instead use list of file 
${PENNCNVPATH}/compile_pfb.pl -listfile ../Samples.txt --snpposfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/pfb/gsa.gr38.pfb


## create gc model file
cd ${PENNCNVPATH}/gc_file/
#gunzip hg38.gc5Base.txt.gz
sort -k 2,2 -k 3,3n hg38.gc5Base.txt > hg38.gc5BaseSorted.txt
cd ${DATADIR}/SNPdata/CNV/PennCNVInput
${PENNCNVPATH}/cal_gc_snp.pl ${PENNCNVPATH}/gc_file/hg38.gc5BaseSorted.txt ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/lib/GSA.gcmodel;


## call CNVs
cd ${DATADIR}/SNPdata/CNV/
mkdir PennCNVOutput

## only run CNVs on samples that pass SNP QC
cd PennCNVInput/
${PENNCNVPATH}/detect_cnv.pl --test -hmm ${PENNCNVPATH}/lib/hh550.hmm -pfb ${PENNCNVPATH}/pfb/gsa.gr38.pfb --listfile ../Samples.txt -log ../PennCNVOutput/SCZ2.log -out ../PennCNVOutput/SCZ_GCModel.rawcnv -gcmodel ${PENNCNVPATH}/lib/GSA.gcmodel;

## generate sample QC report 
${PENNCNVPATH}/filter_cnv.pl --qcsumout PennCNVOutput/SampleQC.txt --qclogfile PennCNVOutput/SCZ2.log PennCNVOutput/SCZ_GCModel.rawcnv




