## Written by Eilis
## Script to perform CNV calling with PennCNV on European samples only
## generally following user guide: https://github.com/WGLab/PennCNV/blob/master/docs/user-guide/input.md
CNVPATH=${DATADIR}/SNPdata/CNV
## assumes executed from script directory
cd ${DATADIR}/
cd SNPdata/Merged

grep "EUR" PredictedPopulations.csv | cut --delim="," -f 2 > EURSamplesTmp.txt
tail -n +2 EURSamplesTmp.txt > EURSamples.txt ## exclude header
rm EURSamplesTmp.txt

## create pfb file
cd ${CNVPATH}/PennCNVInput/

## instead use list of file 
${PENNCNVPATH}/compile_pfb.pl -listfile ${DATADIR}/SNPdata/Merged/EURSamples.txt --snpposfile ${PENNCNVPATH}/lib/gsa.gr38.snpposfile -output ${PENNCNVPATH}/pfb/gsa.gr38.eur.pfb

## call CNVs
cd ${CNVPATH}
cd PennCNVOutput
mkdir EURonly
cd ..

## only run CNVs on samples that pass SNP QC
cd PennCNVInput/
${PENNCNVPATH}/detect_cnv.pl --test -hmm ${PENNCNVPATH}/lib/hh550.hmm -pfb ${PENNCNVPATH}/pfb/gsa.gr38.eur.pfb --listfile ${DATADIR}/SNPdata/Merged/EURSamples.txt -log ../PennCNVOutput/EURonly/SCZ2.log -out ../PennCNVOutput/EURonly/SCZ_GCModel.rawcnv -gcmodel ${PENNCNVPATH}/lib/GSA.gcmodel;

## generate sample QC report 
${PENNCNVPATH}/filter_cnv.pl --qcsumout ../PennCNVOutput/EURonly/SampleQC.txt --qclogfile ../PennCNVOutput/EURonly/SCZ2.log ../PennCNVOutput/EURonly/SCZ_GCModel.rawcnv





