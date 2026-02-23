#!/bin/bash -l
#SBATCH -A project_name   # Replace with your SLURM allocation
#SBATCH -J fastqc_array  
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --cpus-per-task=32
#SBATCH -t 05:00:00
#SBATCH --array=0-80%10
#SBATCH --output=logs/fastqc_%A_%a.out
#SBATCH --mail-user=your_email@domain.com   # Replace with your email
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load fastqc/0.12.1

# USER CONFIGURATION (edit these paths)
RAW_FASTQ_DIR="/path/to/raw_fastq_directory"
SAMPLE_LIST="/path/to/sample_list.txt"               # one sample ID per line
FASTQC_RESULTS_DIR="/path/to/fastqc_output_directory"

mkdir -p "$FASTQC_RESULTS_DIR"

# Retrieve sample ID for this array task
# SLURM array index is 0-based; sed line numbers are 1-based → add +1
ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")

# Exit if no sample was found
if [[ -z "${ACC:-}" ]]; then
  echo "ERROR: No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 1
fi
echo "Running FastQC for sample: $ACC"

# Locate paired FASTQ files
R1=$(find "$RAW_FASTQ_DIR" -type f -name "${ACC}_*_R1_001.fastq.gz" | head -n1)
R2=$(find "$RAW_FASTQ_DIR" -type f -name "${ACC}_*_R2_001.fastq.gz" | head -n1)

if [[ -z "${R1:-}" || -z "${R2:-}" ]]; then
  echo "ERROR: FASTQ files not found for sample $ACC"
  exit 1
fi

# Node-local scratch (PDC_TMP)
TMPDIR="${PDC_TMP:-/tmp}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd "$TMPDIR"

# Copy reads to scratch and run FastQC locally
cp "$R1" "$R2" .
R1_FILE=$(basename "$R1")
R2_FILE=$(basename "$R2")

fastqc -t "${SLURM_CPUS_PER_TASK}" -o "$TMPDIR" "$R1_FILE" "$R2_FILE"

# Copy outputs back
cp -f "$TMPDIR"/*_fastqc.zip "$FASTQC_RESULTS_DIR/" 2>/dev/null || true
cp -f "$TMPDIR"/*_fastqc.html "$FASTQC_RESULTS_DIR/" 2>/dev/null || true

echo "Finished FastQC for sample: $ACC"

