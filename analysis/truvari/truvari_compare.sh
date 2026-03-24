#!/bin/bash
#SBATCH -A naiss2025-22-1574
#SBATCH -p shared
#SBATCH -t 01:00:00
#SBATCH -J truvari_compare
#SBATCH --cpus-per-task=16
#SBATCH -o truvari_compare_%j.out
#SBATCH -e truvari_compare_%j.err
#SBATCH --mail-user=<MY_EMAIL>
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bioinfo-tools
module load bcftools
module load htslib

# Main working directory
WORKDIR="/path/to/project/remapping/merged_vcf/TRUVARI"

# Truvari environment
TRUVARI_ENV="/path/to/truvari_env"
source "${TRUVARI_ENV}/bin/activate"

# Reference FASTA
REF="/path/to/reference/Salix_OZ.fa"

# Input VCFs
LINEAR="${WORKDIR}/salix_linear_PASS_SV50.vcf.gz"
GRAPH="${WORKDIR}/salix_graph_merged_PASS_SV50.vcf.gz"

# Output directory
OUTDIR="${WORKDIR}/truvari_comparison"

echo "=== Input files ==="
echo "Linear VCF: $LINEAR"
echo "Graph VCF:  $GRAPH"
echo "Reference:  $REF"
echo "Workdir:    $WORKDIR"
echo "Outdir:     $OUTDIR"

test -d "$WORKDIR" || { echo "ERROR: Missing WORKDIR: $WORKDIR"; exit 1; }
test -f "$LINEAR" || { echo "ERROR: Missing linear VCF: $LINEAR"; exit 1; }
test -f "${LINEAR}.tbi" || { echo "ERROR: Missing index: ${LINEAR}.tbi"; exit 1; }
test -f "$GRAPH" || { echo "ERROR: Missing graph VCF: $GRAPH"; exit 1; }
test -f "${GRAPH}.tbi" || { echo "ERROR: Missing index: ${GRAPH}.tbi"; exit 1; }
test -f "$REF" || { echo "ERROR: Missing reference: $REF"; exit 1; }

if [[ -d "$OUTDIR" ]]; then
    echo "Removing old Truvari output directory: $OUTDIR"
    rm -rf "$OUTDIR"
fi

echo "=== Truvari help check ==="
truvari bench -h | head -20 || true

echo "=== Running Truvari bench ==="
truvari bench \
    -b "$LINEAR" \
    -c "$GRAPH" \
    -f "$REF" \
    -o "$OUTDIR" \
    -p 0.7 \
    -P 0.7 \
    -s 50 \
    --passonly

echo "=== Truvari finished ==="

echo "=== Summary ==="
cat "$OUTDIR/summary.txt"

echo "=== Output files ==="
ls -lh "$OUTDIR"





#!/bin/bash
#SBATCH -A naiss2025-22-1574
#SBATCH -p shared
#SBATCH -t 00:01:00
#SBATCH -J truvari_compare
#SBATCH --cpus-per-task=16
#SBATCH -o truvari_compare_%j.out
#SBATCH -e truvari_compare_%j.err
#SBATCH --mail-user=<MY_EMAIL>
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bioinfo-tools
module load bcftools
module load htslib

# Truvari environment
TRUVARI_ENV="/path/to/truvari_env"
source "${TRUVARI_ENV}/bin/activate"

# Main working directory
WORKDIR="/path/to/project/remapping/merged_vcf/TRUVARI"

# Truvari environment
TRUVARI_ENV="${WORKDIR}"
source "${TRUVARI_ENV}/bin/activate"

# Reference FASTA
REF="/path/to/reference/Salix_OZ.fa"

# Input VCFs
LINEAR="${WORKDIR}/salix_linear_PASS_SV50.vcf.gz"
GRAPH="${WORKDIR}/salix_graph_merged_PASS_SV50.vcf.gz"

# Output directory
OUTDIR="${WORKDIR}/truvari_comparison"

# SANITY CHECKS
echo "=== Input files ==="
echo "Linear VCF: $LINEAR"
echo "Graph VCF:  $GRAPH"
echo "Reference:  $REF"
echo "Workdir:    $WORKDIR"
echo "Outdir:     $OUTDIR"

test -d "$WORKDIR" || { echo "ERROR: Missing WORKDIR: $WORKDIR"; exit 1; }
test -f "$LINEAR" || { echo "ERROR: Missing linear VCF: $LINEAR"; exit 1; }
test -f "${LINEAR}.tbi" || { echo "ERROR: Missing index: ${LINEAR}.tbi"; exit 1; }
test -f "$GRAPH" || { echo "ERROR: Missing graph VCF: $GRAPH"; exit 1; }
test -f "${GRAPH}.tbi" || { echo "ERROR: Missing index: ${GRAPH}.tbi"; exit 1; }
test -f "$REF" || { echo "ERROR: Missing reference: $REF"; exit 1; }

# PREPARE OUTPUT DIRECTORY
if [[ -d "$OUTDIR" ]]; then
    echo "Removing old Truvari output directory: $OUTDIR"
    rm -rf "$OUTDIR"
fi

# CHECK TRUVARI
echo "=== Truvari help check ==="
truvari bench -h | head -20 || true

# RUN TRUVARI
echo "=== Running Truvari bench ==="
truvari bench \
    -b "$LINEAR" \
    -c "$GRAPH" \
    -f "$REF" \
    -o "$OUTDIR" \
    -p 0.7 \
    -P 0.7 \
    -s 50 \
    --passonly

echo "=== Truvari finished ==="

echo "=== Summary ==="
cat "$OUTDIR/summary.txt"

echo "=== Output files ==="
ls -lh "$OUTDIR"
