#!/bin/bash
#SBATCH -A naiss2025-22-1574
#SBATCH -p shared
#SBATCH -t 00:10:00
#SBATCH -J vg_merge_vcf
#SBATCH --cpus-per-task=32
#SBATCH -o ${PROJECT_DIR}/logs/%x_%j.out
#SBATCH -e ${PROJECT_DIR}/logs/%x_%j.err
#SBATCH --mail-user=<MY_EMAIL>        # replace
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bioinfo-tools
module load bcftools
module load htslib

# Main project directory 
PROJECT_DIR="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/remapping"

# Input/output directories
VCFDIR="${PROJECT_DIR}/genotyped_vcf"
OUTDIR="${PROJECT_DIR}/merged_vcf"
LOGDIR="${PROJECT_DIR}/logs"

# Output file
OUTVCF="${OUTDIR}/salix_graph_samples.merged.vcf.gz"

# Threads for bcftools
THREADS=4

mkdir -p "$OUTDIR" "$LOGDIR"

FILELIST="${OUTDIR}/vcf_list.txt"

echo "== Sanity checks =="
test -d "$VCFDIR" || { echo "ERROR: Missing VCFDIR: $VCFDIR"; exit 1; }

# COLLECT VCF FILES
echo "== Collecting input VCFs =="
find "$VCFDIR" -name "*.genotyped.vcf.gz" | sort > "$FILELIST"

NVCF=$(wc -l < "$FILELIST")
echo "Found $NVCF compressed per-sample VCF files"

if [[ "$NVCF" -eq 0 ]]; then
    echo "ERROR: no *.genotyped.vcf.gz files found in $VCFDIR"
    exit 1
fi

# ENSURE INDEXES
echo "== Checking indexes =="
while read -r f; do
    [[ -f "${f}.tbi" ]] || tabix -f -p vcf "$f"
done < "$FILELIST"

# MERGE VCFs
echo "== Merging VCFs with bcftools =="
bcftools merge \
    --threads "$THREADS" \
    -l "$FILELIST" \
    -Oz \
    -o "$OUTVCF"

# INDEX MERGED VCF
echo "== Index merged VCF =="
tabix -f -p vcf "$OUTVCF"

echo "== Done =="
ls -lh "$OUTVCF" "${OUTVCF}.tbi"

echo "Number of samples in merged VCF:"
bcftools query -l "$OUTVCF" | wc -l
