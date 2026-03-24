#!/bin/bash
#SBATCH -A naiss2025-22-1574                  
#SBATCH -J liftoff
#SBATCH -t 12:00:00
#SBATCH -p shared
#SBATCH --cpus-per-task=24
#SBATCH -o liftoff_%j.out
#SBATCH -e liftoff_%j.err
#SBATCH --mail-user=<MY_EMAIL> 
#SBATCH --mail-type=END,FAIL

set -euo pipefail

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# Project directory 
PROJECT_DIR="/path/to/reference"

# Input files
TARGET_GENOME="${PROJECT_DIR}/Salix_viminalis.fasta"
REFERENCE_GENOME="${PROJECT_DIR}/osier_cornell.fa"
REFERENCE_GFF="${PROJECT_DIR}/osier_cornell.gff"

# Output directory
OUTDIR="${PROJECT_DIR}/liftoff_out"

# Conda environment
CONDA_BASE="/path/to/miniconda3"
CONDA_ENV="liftoff_env"

# Threads
THREADS=24

# SETUP
mkdir -p "$OUTDIR"

echo "=== Input files ==="
echo "Target genome:    $TARGET_GENOME"
echo "Reference genome: $REFERENCE_GENOME"
echo "Reference GFF:    $REFERENCE_GFF"
echo "Output dir:       $OUTDIR"

test -f "$TARGET_GENOME" || { echo "ERROR: Missing target genome"; exit 1; }
test -f "$REFERENCE_GENOME" || { echo "ERROR: Missing reference genome"; exit 1; }
test -f "$REFERENCE_GFF" || { echo "ERROR: Missing reference GFF"; exit 1; }

# RUN LIFTOFF
echo "=== Running Liftoff ==="

liftoff "$TARGET_GENOME" "$REFERENCE_GENOME" \
  -g "$REFERENCE_GFF" \
  -o "$OUTDIR/lifted.gff3" \
  -dir "$OUTDIR" \
  -p "$THREADS" \
  -copies \
  -polish \
  -sc 0.9

echo "=== Liftoff completed ==="

ls -lh "$OUTDIR"
