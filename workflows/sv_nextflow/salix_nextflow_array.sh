#!/bin/bash
#SBATCH -A naiss2025-22-1574
#SBATCH -J sv_nf
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -t 5-00:00:00
#SBATCH --output=/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/logs/sv_pipeline_all_bams%j.out
#SBATCH --mail-user=saida.sharifova.1623@student.uu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module purge
module load bcftools/1.20
module load bedtools/2.31.0
module load nextflow
module load delly/1.1.6
module load manta/1.6.0
module load samtools/1.20
module load cray-python/3.11.7
module load dysgu/1.8.7

export PATH="$HOME/.local/bin:$PATH"

export NXF_HOME="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Nextflow_outputs/.nextflow"
export NXF_WORK="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Nextflow_outputs/Nextflow_work"
mkdir -p "$NXF_HOME" "$NXF_WORK"

export NXF_OPTS="-Xms1g -Xmx4g"

nextflow run sv_calling.nf \
  -c salix_nextflow_array.config \
  -work-dir "$NXF_WORK" \
  --reference '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Reference/Salix_viminalis.fasta' \
  --bam_dir '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Mapping_output/bam_salix_all/' \
  --chr_list '/cfs/klemming/projects/supr/naiss2025-23-666/Saida/chr_list.txt' \
  -with-trace trace_array.txt \
  -resume
