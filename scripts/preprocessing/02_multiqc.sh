#!/bin/bash -l
#SBATCH -A project_name    # Replace with your SLURM allocation
#SBATCH -J multiqc_fastqc
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -t 05:00:00
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --mail-user=your_email@domain.com   # Replace with your email
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load multiqc/1.30

# USER CONFIGURATION    (edit these paths)
FASTQC_RESULTS_DIR="/path/to/fastqc_output_directory"     # contains *_fastqc.zip/html
MULTIQC_OUT_DIR="/path/to/multiqc_output_directory"

mkdir -p "$MULTIQC_OUT_DIR"

TMPDIR="${PDC_TMP:-/tmp}/multiqc_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# Copy inputs to scratch
cp -r "$FASTQC_RESULTS_DIR"/. "$TMPDIR/"

echo "Running MultiQC in node-local TMPDIR..."
multiqc "$TMPDIR" -o "$TMPDIR"

# Copy outputs back
cp -f "$TMPDIR/multiqc_report.html" "$MULTIQC_OUT_DIR/"
cp -r "$TMPDIR/multiqc_data" "$MULTIQC_OUT_DIR/" 2>/dev/null || true

echo "MultiQC report generated:"
echo "$MULTIQC_OUT_DIR/multiqc_report.html"

