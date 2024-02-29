#!/bin/bash

# ============================================================================ #
# Uses LD score regression software (ldsc) with cusutom annotation derived from
# genomics data (like 1000 genomes) to estimate proportion of heritability
# attributitable to SNPs in these regions 
# Example reference files can be downloaded from: 
# https://zenodo.org/records/10515792

# NB: although annotations are on hg38, GWAS files don't have location 
# information so are matched by rsID hence no issue with genome builds 
# for GWAS traits
# ============================================================================ #

## Compute LD scores with custom annotation file.
for chr in {1..22};
do
  python \
	"${LD_SOFTWARE_DIR}/ldsc.py" \
	--l2 \
	--bfile      "${LD_REFERENCE_DIR}/plink_files/${REFERENCE_PREFIX}.${chr}" \
  --ld-wind-cm 1 \
  --annot      "${LD_ANNOTATION_DIR}/${ANNOTATION_PREFIX}.${chr}.annot" \
  --out        "${LD_ANNOTATION_DIR}/${ANNOTATION_PREFIX}/${ANNOTATION_PREFIX}.${chr}" \
  --print-snps "${SNP_LISTS_DIR}/${SNP_LIST_PREFIX}.${chr}.snp"
done

## estimate partioned heritability for a range of GWAS traits
gwastraits=$(ls "${LD_GWAS_TRAITS_DIR}/*${GWAS_PATTERN}*.gz")

for filename in "${gwastraits[@]}"; 
do
  outfile=$(basename "${filename}" .sumstats.gz)

	python \
	"${LD_SOFTWARE_DIR}/ldsc.py" \
	--h2          "${filename}" \
	--ref-ld-chr  "${LD_ANNOTATION_DIR}/${ANNOTATION_PREFIX}/${ANNOTATION_PREFIX}." \
	--frqfile-chr "${LD_REFERENCE_DIR}/frq_files/${REFERENCE_PREFIX}." \
	--w-ld-chr    "${LD_REFERENCE_DIR}/weights/${WEIGHTS_PREFIX}." \
	--out         "${OUTPUTS_DIR}/${ANNOTATION_PREFIX}/${outfile}" \
	--overlap-annot \
	--print-coefficients
done
