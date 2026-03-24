#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Filter SV VCF by number of called samples and extract inputs for downstream PCA / R analysis

# Input:
#   final_overlap.vcf.gz
#
# Output:
#   header.vcf
#   body.vcf
#   sv_called10.vcf.gz
#   sv_called10.vcf.gz.tbi
#   sv_GT.tsv
#   samples.txt
#
# Description:
#   1. Extract VCF header
#   2. Keep only variants called in at least 10 samples (genotype field not equal to ./.)
#   3. Rebuild and compress filtered VCF
#   4. Index filtered VCF
#   5. Extract genotype table for R
#   6. Extract sample names
# ============================================================

# Input file
VCF="final_overlap.vcf.gz"

# Output files

HEADER="header.vcf"
BODY="body.vcf"
FILTERED_VCF="sv_called10.vcf"
GT_TABLE="sv_GT.tsv"
SAMPLE_LIST="samples.txt"

# Extract VCF header
bcftools view -h "${VCF}" > "${HEADER}"

# Filter variants
# Keep only records with >=10 called samples
# - Sample genotype columns start at column 10
# - We count a sample as "called" if GT != ./.
# - Assumes GT is the first field in each sample column

bcftools view -H "${VCF}" | \
awk -F'\t' 'BEGIN{OFS="\t"} {
  called = 0
  for (i = 10; i <= NF; i++) {
    split($i, a, ":")
    if (a[1] != "./.") called++
  }
  if (called >= 10) print
}' > "${BODY}"

# Rebuild filtered VCF
cat "${HEADER}" "${BODY}" > "${FILTERED_VCF}"

# Compress and index filtered VCF
bgzip -f "${FILTERED_VCF}"
tabix -p vcf "${FILTERED_VCF}.gz"

# Extract genotype table for R
# Output columns:
#   CHROM  POS  ID  sample1_GT  sample2_GT ...

bcftools query \
  -f '%CHROM\t%POS\t%ID[\t%GT]\n' \
  "${FILTERED_VCF}.gz" > "${GT_TABLE}"

# Extract sample names
bcftools query -l "${FILTERED_VCF}.gz" > "${SAMPLE_LIST}"

# Done
echo "Filtering complete."
echo "Filtered VCF: ${FILTERED_VCF}.gz"
echo "Genotype table: ${GT_TABLE}"
echo "Sample list: ${SAMPLE_LIST}"
