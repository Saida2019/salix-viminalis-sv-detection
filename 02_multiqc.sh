#!/bin/bash -l
#SBATCH -A project_name
#SBATCH -J multiqc_fastqc
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -t 05:00:00
#SBATCH --output=logs/multiqc_%j.out
#SBATCH --mail-user=your_email@domain.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load multiqc/1.30

# USER CONFIGURATION (edit these paths)
FASTQC_RESULTS_DIR="/path/to/fastqc_output_directory"     # contains *_fastqc.zip/html
MULTIQC_OUT_DIR="/path/to/multiqc_output_directory"

mkdir -p "$MULTIQC_OUT_DIR"

TMPDIR="${PDC_TMP:-/tmp}/multiqc_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# Copy inputs to scratch (helps performance)
cp -r "$FASTQC_RESULTS_DIR"/. "$TMPDIR/"

echo "Running MultiQC in node-local TMPDIR..."
multiqc "$TMPDIR" -o "$TMPDIR"

if [[ ! -f "$TMPDIR/multiqc_report.html" ]]; then
  echo "ERROR: multiqc_report.html not found — MultiQC may have failed."
  exit 1
fi

# Copy report + data folder back
cp -f "$TMPDIR/multiqc_report.html" "$MULTIQC_OUT_DIR/"
if [[ -d "$TMPDIR/multiqc_data" ]]; then
  rm -rf "$MULTIQC_OUT_DIR/multiqc_data" 2>/dev/null || true
  cp -r "$TMPDIR/multiqc_data" "$MULTIQC_OUT_DIR/"
fi

echo "MultiQC report generated:"
echo "$MULTIQC_OUT_DIR/multiqc_report.html"

