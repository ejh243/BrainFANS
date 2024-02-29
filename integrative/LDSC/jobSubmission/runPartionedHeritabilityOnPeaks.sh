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

module purge
module load Anaconda3/2020.02

source "${CONDA_SHELL}/profile.d/conda.sh" || \
{ echo "profile.d/conda.sh does not exist in specified location: \
[\${CONDA_SHELL} - ${CONDA_SHELL}]"; exit 1; }
conda activate "${LDSC_CONDA_ENVIRONMENT}"

## Compute LD scores with custom annotation file.
for chr in {1..22}; do
  python \
	"${LD_SOFTWARE_DIR}/ldsc.py" \
	--l2 \
	--bfile      "${LD_REFERENCE_DIR}/plink_files/${REFERENCE_PREFIX}.${chr}" \
  --ld-wind-cm 1 \
  --annot      "${LD_ANNOTATION_DIR}/${ANNOTATION_PREFIX}.${chr}.annot" \
  --out        "${OUTPUTS_DIR}/${ANNOTATION_PREFIX}/${ANNOTATION_PREFIX}.${chr}" \
  --print-snps "${SNP_LISTS_DIR}/${SNP_LIST_PREFIX}.${chr}.snp"
done

## estimate partioned heritability for a given selection of GWAS traits
gwas_traits=$(find "${LD_GWAS_TRAITS_DIR}" -name "*${GWAS_PATTERN}*.gz")

for file_name in "${gwas_traits[@]}"; do
  output_file=$(basename "${file_name}" .sumstats.gz)

	python \
	"${LD_SOFTWARE_DIR}/ldsc.py" \
	--h2          "${file_name}" \
	--ref-ld-chr  "${OUTPUTS_DIR}/${ANNOTATION_PREFIX}/${ANNOTATION_PREFIX}." \
	--frqfile-chr "${LD_REFERENCE_DIR}/frq_files/${REFERENCE_PREFIX}." \
	--w-ld-chr    "${LD_REFERENCE_DIR}/weights/${WEIGHTS_PREFIX}." \
	--out         "${OUTPUTS_DIR}/${ANNOTATION_PREFIX}/heritability/${output_file}" \
	--overlap-annot \
	--print-coefficients
done
