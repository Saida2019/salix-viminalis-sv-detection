#!/bin/bash -l
#SBATCH -A project_name
#SBATCH -J addRG_array_picard
#SBATCH -p shared
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --array=0-80%5
#SBATCH --output=logs/rg_%A_%a.out
#SBATCH --mail-user=your_email@domain.com
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load java/17.0.4
module load picard/3.3.0

# USER CONFIGURATION (edit these paths)
BAM_IN_DIR="/path/to/input/bams"                    
BAM_RG_DIR="/path/to/output/bams_with_rg"
SAMPLE_TABLE="/path/to/samples_bam_table.txt" 

# Expected format of SAMPLE_TABLE:
# Header + one BAM per line
# Columns: bam_filename  SM  LB  PL  PU  PM  PI
# Note: Only the first column (bam_filename) is used in this script.

mkdir -p "$BAM_RG_DIR"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 2))p" "$SAMPLE_TABLE")
if [[ -z "${LINE:-}" ]]; then
  echo "ERROR: No line found in SAMPLE_TABLE for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

read -r BAM SM LB PL PU PM PI <<< "$LINE"

IN_BAM="$BAM_IN_DIR/$BAM"
if [[ ! -f "$IN_BAM" ]]; then
  echo "ERROR: Input BAM not found: $IN_BAM"
  exit 1
fi

OUT_BAM="$BAM_RG_DIR/${BAM%.bam}.rg.bam"

echo "Adding read groups to: $IN_BAM (SM=$SM)"

java -jar "$EBROOTPICARD/picard.jar" AddOrReplaceReadGroups \
  I="$IN_BAM" \
  O="$OUT_BAM" \
  RGID="$SM" \
  RGLB="$LB" \
  RGPL="$PL" \
  RGPU="$PU" \
  RGSM="$SM" \
  RGPM="$PM" \
  RGPI="$PI" \
  CREATE_INDEX=true

echo "Finished RG: $OUT_BAM"

