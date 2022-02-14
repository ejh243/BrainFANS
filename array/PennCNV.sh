## Written by Eilis
## Script to perform ${CNVPATH} calling with Penn${CNVPATH}
## generally following user guide: https://github.com/WGLab/Penn${CNVPATH}/blob/master/docs/user-guide/input.md

## assumes executed from script directory
cd ${DATADIR}/
cd SNPdata/Merged
CNVPATH=${DATADIR}/SNPdata/CNV

## write list of samples that passed QC for CNV calling
cut -f 3,4 --delimiter=" " UpdateIDs.txt > ${CNVPATH}/ID_Map.txt
cut -f 2 --delimiter=" " SCZ2_QCd.fam > ${CNVPATH}/Samples.txt

## create input files from final report
## replace array ids with sample IDs nb change "/" within names to "_"
${PENNCNVPATH}/split_illumina_report.pl -prefix ${CNVPATH}/PennCNVInput/ -revised_file ${CNVPATH}/ID_Map.txt FANsBrainsFinalReport.txt
## need to replace column headers
FILES=${CNVPATH}/PennCNVInput/*
for f in $FILES
do
	sed -i '1s/.*/Name\tLog R Ratio\tB Allele Freq/' ${f}
done

##create snp position file excluding markers without a position
awk -v OFS='\t' '{if ($1 != 0) print $2, $1,$4}' SCZ_gr38_binary.bim > ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp
echo -e "Name\tChr\tPos" | cat - ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp > ${PENNCNVPATH}/lib/gsa.gr38.snpposfile
rm ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp


## create pfb files
cd ${PENNCNVPATH}
mkdir -p ./pfb
cd ${CNVPATH}/PennCNVInput/

## instead use list of file 
${PENNCNVPATH}/compile_pfb.pl -listfile ${CNVPATH}/Samples.txt --snpposfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/pfb/gsa.gr38.pfb


## create gc model file
cd ${PENNCNVPATH}/gc_file/
#gunzip hg38.gc5Base.txt.gz
sort -k 2,2 -k 3,3n hg38.gc5Base.txt > hg38.gc5BaseSorted.txt
cd ${CNVPATH}/PennCNVInput
${PENNCNVPATH}/cal_gc_snp.pl ${PENNCNVPATH}/gc_file/hg38.gc5BaseSorted.txt ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/lib/GSA.gcmodel;


## call ${CNVPATH}s
cd ${CNVPATH}/
mkdir -p PennCNVOutput

## only run ${CNVPATH}s on samples that pass SNP QC
cd PennCNVInput/
${PENNCNVPATH}/detect_cnv.pl --test -hmm ${PENNCNVPATH}/lib/hh550.hmm -pfb ${PENNCNVPATH}/pfb/gsa.gr38.pfb --listfile ${CNVPATH}/Samples.txt -log ${CNVPATH}/PennCNVOutput/Merged.log -out ${CNVPATH}/PennCNVOutput/SCZ_GCModel.rawcnv -gcmodel ${PENNCNVPATH}/lib/GSA.gcmodel;

## generate sample QC report 
${PENNCNVPATH}/filter_cnv.pl --qcsumout PennCNVOutput/SampleQC.txt --qclogfile PennCNVOutput/Merged.log PennCNVOutput/SCZ_GCModel.rawcnv




