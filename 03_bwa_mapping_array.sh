#!/bin/bash -l
#SBATCH -A project_name
#SBATCH -J bwa_mapping_array
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH --array=0-80%5                 
#SBATCH --output=logs/bwa_%A_%a.out
#SBATCH --mail-user=your_email@domain.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bwa/0.7.18
module load samtools/1.20

# USER CONFIGURATION (edit these paths)
RAW_FASTQ_DIR="/path/to/raw_fastq_directory"
SAMPLE_LIST="/path/to/sample_list.txt"               # one sample ID per line
REFERENCE_FASTA="/path/to/reference/genome.fasta"
BWA_RESULTS_DIR="/path/to/output/bwa_results_directory"

mkdir -p "$BWA_RESULTS_DIR"

ACC=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")
if [[ -z "${ACC:-}" ]]; then
  echo "ERROR: No sample found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  exit 1
fi
echo "Processing sample: $ACC"

R1=$(find "$RAW_FASTQ_DIR" -type f -name "${ACC}_*_R1_001.fastq.gz" | head -n1)
R2=$(find "$RAW_FASTQ_DIR" -type f -name "${ACC}_*_R2_001.fastq.gz" | head -n1)
if [[ -z "${R1:-}" || -z "${R2:-}" ]]; then
  echo "ERROR: FASTQ files missing for $ACC"
  exit 1
fi

TMPDIR="${PDC_TMP:-/tmp}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
cd "$TMPDIR"

cp "$R1" "$R2" "$TMPDIR"
R1_FILE=$(basename "$R1")
R2_FILE=$(basename "$R2")

BAM_OUT="${ACC}.sorted.bam"

bwa mem -M -t "${SLURM_CPUS_PER_TASK}" "$REFERENCE_FASTA" "$R1_FILE" "$R2_FILE" | \
  samtools view -b -@ "${SLURM_CPUS_PER_TASK}" - | \
  samtools sort -@ "${SLURM_CPUS_PER_TASK}" -m 1G -T "$TMPDIR/sort_tmp" -o "$BAM_OUT" -

samtools index "$BAM_OUT"

mv "$BAM_OUT" "${BAM_OUT}.bai" "$BWA_RESULTS_DIR/"

echo "Finished mapping: $ACC"


