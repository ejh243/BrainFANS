## Written by Eilis
## Script to perform CNV calling with PennCNV
## generally following user guide: https://github.com/WGLab/PennCNV/blob/master/docs/user-guide/input.md

## assumes executed from script directory
mkdir -p ${CNVDIR}/PennCNVInput/

## create input files from final report
${PENNCNVPATH}/split_illumina_report.pl -prefix ${CNVDIR}/PennCNVInput/ ${RAWDATADIR}/${FILEPREFIX}FinalReport.txt
## replace array ids with sample IDs nb change "/" within names to "_"
while read p; do
  old=$(echo $p | cut --delim=" " -f 1)
  new=$(echo $p | cut --delim=" " -f 2)
  if [ -f ${CNVDIR}/PennCNVInput/$old ];
  then
    echo move $old to $new
	mv ${CNVDIR}/PennCNVInput/$old ${CNVDIR}/PennCNVInput/$new
  else
    echo $old already changed
  fi
done < ${CNVDIR}/ID_Map.txt

## need to replace column headers
FILES=${CNVDIR}/PennCNVInput/*
for f in $FILES
do
	sed -i '1s/.*/Name\tLog R Ratio\tB Allele Freq/' ${f}
done

##create snp position file excluding markers without a position
awk -v OFS='\t' '{if ($1 != 0) print $2, $1,$4}' ${RAWDATADIR}/${FILEPREFIX}.bim > ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp
echo -e "Name\tChr\tPos" | cat - ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp > ${PENNCNVPATH}/lib/gsa.gr38.snpposfile
rm ${PENNCNVPATH}/lib/gsa.gr38.snpposfile.tmp


## create pfb files
cd ${PENNCNVPATH}
mkdir -p ./pfb
cd ${CNVDIR}/PennCNVInput/

## instead use list of file 
${PENNCNVPATH}/compile_pfb.pl -listfile ${CNVDIR}/Samples.txt --snpposfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/pfb/gsa.gr38.pfb


## create gc model file
cd ${PENNCNVPATH}/gc_file/
#gunzip hg38.gc5Base.txt.gz
sort -k 2,2 -k 3,3n hg38.gc5Base.txt > hg38.gc5BaseSorted.txt
cd ${CNVDIR}/PennCNVInput
${PENNCNVPATH}/cal_gc_snp.pl ${PENNCNVPATH}/gc_file/hg38.gc5BaseSorted.txt ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/lib/GSA.gcmodel;


## call cnvs
cd ${CNVDIR}/
mkdir -p PennCNVOutput

## only run calling on samples that pass SNP QC
cd PennCNVInput/
${PENNCNVPATH}/detect_cnv.pl --test -hmm ${PENNCNVPATH}/lib/hh550.hmm -pfb ${PENNCNVPATH}/pfb/gsa.gr38.pfb --listfile ${CNVDIR}/Samples.txt -log ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_Merged.log -out ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_GCModel.rawcnv -gcmodel ${PENNCNVPATH}/lib/GSA.gcmodel;

## generate sample QC report 
${PENNCNVPATH}/filter_cnv.pl --qcsumout ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_SampleQC.txt --qclogfile ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_Merged.log ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_GCModel.rawcnv




