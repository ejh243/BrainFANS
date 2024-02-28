#!/bin/bash

## use LD score regression software with cusutom annotaion derived from genomics data to estimate proportion of heritability attributitable to SNPs in these regions 
## uses reference files downloaded from https://data.broadinstitute.org/alkesgroup/LDSCORE/
## NB although annotations are on hg38, GWAS files don't have location information so are matched by rsID hence no issue with genome builds for GWAS traits

## Compute LD scores with custom annot file.
for chr in {1..22};
do
  python ${LDPATH}/ldsc.py --l2 --bfile ${LDREFPATH}grch38/plink_files/1000G.EUR.hg38.${chr}\
  --ld-wind-cm 1\
  --annot ${LDANNOPATH}NeuralCellRegulatoryPeaks.${chr}.annot\
  --out ${LDANNOPATH}NeuralCellRegulatoryPeaks.${chr}\
  --print-snps ${LDREFPATH}/hapmap3_snps/hm.${chr}.snp
done

## estimate partioned heritability for a range of GWAS traits
gwastraits=($ls ${LDREFPATH}/gwas_traits/Anorexia*.gz)

for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDANNOPATH}NeuralCellRegulatoryPeaks. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile \
	--print-coefficients
done

## just limit to Neural peak sets to see if results change

for chr in {1..22};
do
  cut -f 1-13 -d' ' ${LDANNOPATH}NeuralCellRegulatoryPeaks.${chr}.annot > ${LDANNOPATH}NeuralCellRegulatoryPeaksonly.${chr}.annot
  python ${LDPATH}/ldsc.py --l2 --bfile ${LDREFPATH}grch38/plink_files/1000G.EUR.hg38.${chr}\
  --ld-wind-cm 1\
  --annot ${LDANNOPATH}NeuralCellRegulatoryPeaksonly.${chr}.annot\
  --out ${LDANNOPATH}\NeuralCellRegulatoryPeaksonly.${chr}\
  --print-snps ${LDREFPATH}/hapmap3_snps/hm.${chr}.snp
done

for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDANNOPATH}\NeuralCellRegulatoryPeaksonly. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile.NeuralCellRegulatoryPeaksOnly \
	--print-coefficients
done



for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDANNOPATH}/NeuralCellRegulatoryPeaksonly.,${LDREFPATH}/grch38/baselineLD_v2.2/Rerun/baselineLD.Rerun. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile.NeuralCellRegulatoryPeaksVsBaselinev2.2 \
	--print-coefficients
done

