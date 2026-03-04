#!/bin/bash
#SBATCH -A naiss2025-22-1574
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH -J salix_vg_pack_call
#SBATCH --cpus-per-task=96
#SBATCH --array=1-80%3
#SBATCH -o /cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/remapping/logs/%x_%A_%a.out
#SBATCH -e /cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/remapping/logs/%x_%A_%a.err
#SBATCH --mail-user="my_email" #change
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bioinfo-tools
module load vg/1.48.0
module load htslib

WORK="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/remapping"

# Indexes
XG="$WORK/salix.sv.xg"

# SV VCF (the same used to build the graph)
SV_VCF="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/graph_pangenome/vcf.noDUPTANDEM_keepINV.vcf.gz"

# Sample sheet (filtered 80-sample list)
SAMPLES_TSV="$WORK/samples_80.tsv"

# Inputs (GAMs) from mapping step
GAMDIR="$WORK/mapping_gam"

# Outputs
PACKDIR="$WORK/packed"
VCFDIR="$WORK/genotyped_vcf"
mkdir -p "$WORK/logs" "$PACKDIR" "$VCFDIR"

THREADS=16
MQ=5

echo "== Sanity checks =="
test -f "$XG" || { echo "Missing XG: $XG"; exit 1; }
test -f "$SV_VCF" || { echo "Missing SV_VCF: $SV_VCF"; exit 1; }
test -f "${SV_VCF}.tbi" || { echo "Missing VCF index: ${SV_VCF}.tbi"; exit 1; }
test -f "$SAMPLES_TSV" || { echo "Missing sample list: $SAMPLES_TSV"; exit 1; }

# Grab the array line
LINE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TSV" || true)"
if [[ -z "${LINE}" ]]; then
  echo "No line found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in $SAMPLES_TSV"
  exit 1
fi

# Sample ID is first column; strip CR in case file has Windows line endings
SAMPLE_ID="$(echo "$LINE" | cut -f1 | tr -d '\r')"
test -n "$SAMPLE_ID" || { echo "Empty SAMPLE_ID on line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

GAMGZ="$GAMDIR/${SAMPLE_ID}.vs_graph.gam.gz"
PACK="$PACKDIR/${SAMPLE_ID}.vs_graph.pack"
OUTVCF="$VCFDIR/${SAMPLE_ID}.genotyped.vcf"

echo "== Task info =="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE_ID}"
echo "Threads: ${THREADS}"
echo "GAMGZ: ${GAMGZ}"
echo "PACK: ${PACK}"
echo "OUTVCF: ${OUTVCF}"

test -f "$GAMGZ" || { echo "Missing GAM.GZ: $GAMGZ"; exit 1; }

echo "== Step A: vg pack (read support) =="
# -Q filters low-MQ alignments
zcat "$GAMGZ" | vg pack -x "$XG" -g - -o "$PACK" -Q "$MQ" -t "$THREADS"

echo "== Step B: vg call (genotype known SV sites) =="
# -v uses the original SV VCF (sites) for genotyping
# -s sets sample name in output VCF
vg call "$XG" -k "$PACK" -s "$SAMPLE_ID" -v "$SV_VCF" -t "$THREADS" > "$OUTVCF"

# bgzip + tabix for downstream tools
if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
  bgzip -f "$OUTVCF"
  tabix -f -p vcf "${OUTVCF}.gz"
  echo "Compressed+indexed: ${OUTVCF}.gz"
else
  echo "NOTE: bgzip/tabix not found in PATH; left as plain VCF."
fi

echo "== Done =="
ls -lh "$PACK" "${OUTVCF}.gz" 2>/dev/null || ls -lh "$PACK" "$OUTVCF"




