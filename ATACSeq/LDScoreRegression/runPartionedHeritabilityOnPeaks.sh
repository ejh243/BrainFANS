## use LD score regression software with cusutom annotaion derived from ATAQ-peaks to estimate proportion of heritability attributitable to SNPs in these regions 
## uses reference files downloaded from https://data.broadinstitute.org/alkesgroup/LDSCORE/


## Compute LD scores with custom annot file.
for chr in {1..22};
do
  python ${LDPATH}/ldsc.py --l2 --bfile ${LDREFPATH}grch38/plink_files/1000G.EUR.hg38.${chr}\
  --ld-wind-cm 1\
  --annot ${LDANNOPATH}ATACPeaks.${chr}.annot\
  --out ${LDANNOPATH}\ATACPeaks.${chr}\
  --print-snps ${LDREFPATH}/hapmap3_snps/hm.${chr}.snp
done

## estimate partioned heritability for a range of GWAS traits
gwastraits=($ls ${LDREFPATH}/gwas_traits/*.gz)

for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDANNOPATH}\ATACPeaks. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile \
	--print-coefficients
done

## just limit to ATAC peak sets to see if results change

for chr in {1..22};
do
  cut -f 1-8 -d' ' ${LDANNOPATH}ATACPeaks.${chr}.annot > ${LDANNOPATH}ATACPeaksOnly.${chr}.annot
  python ${LDPATH}/ldsc.py --l2 --bfile ${LDREFPATH}grch38/plink_files/1000G.EUR.hg38.${chr}\
  --ld-wind-cm 1\
  --annot ${LDANNOPATH}ATACPeaksOnly.${chr}.annot\
  --out ${LDANNOPATH}\ATACPeaksOnly.${chr}\
  --print-snps ${LDREFPATH}/hapmap3_snps/hm.${chr}.snp
done

for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDANNOPATH}\ATACPeaksOnly. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile.ATACPeaksOnly \
	--print-coefficients
done

## run full baseline model for comparision

for chr in {1..22};
do
  python ${LDPATH}/ldsc.py --l2 --bfile ${LDREFPATH}grch38/plink_files/1000G.EUR.hg38.${chr}\
  --ld-wind-cm 1\
  --annot ${LDREFPATH}/grch38/baselineLD_v2.2/baselineLD.${chr}.annot.gz\
  --out ${LDREFPATH}/grch38/baselineLD_v2.2/baselineLD.Rerun.${chr}\
  --print-snps ${LDREFPATH}/hapmap3_snps/hm.${chr}.snp
done

for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDREFPATH}/grch38/baselineLD_v2.2/Rerun/baselineLD.Rerun. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile.Baselinev2.2 \
	--print-coefficients
done

for filename in ${gwastraits[@]}; 
do
    outfile=$(basename $filename .sumstats.gz)
	python ${LDPATH}/ldsc.py --h2 $filename \
	--ref-ld-chr ${LDANNOPATH}/ATACPeaksOnly.,${LDREFPATH}/grch38/baselineLD_v2.2/Rerun/baselineLD.Rerun. \
	--overlap-annot \
	--frqfile-chr ${LDREFPATH}1000G_Phase3_frq/1000G.EUR.QC. \
	--w-ld-chr  ${LDREFPATH}grch38/weights/weights.hm3_noMHC. --out ${OUTPATH}/$outfile.ATACPeaksVsBaselinev2.2 \
	--print-coefficients
done

