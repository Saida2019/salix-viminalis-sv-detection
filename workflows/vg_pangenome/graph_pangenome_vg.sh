#!/bin/bash
#SBATCH -A naiss2025-22-1574
#SBATCH -p shared
#SBATCH -t 96:00:00
#SBATCH -J salix_vg_build
#SBATCH --cpus-per-task=128
#SBATCH -o /cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/graph_pangenome/05_logs/%x_%j.out
#SBATCH -e /cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/graph_pangenome/05_logs/%x_%j.err
#SBATCH --mail-user="my_email" #change
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load bioinfo-tools
module load vg/1.48.0

REF="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Reference/OZ_only/Salix_OZ.fa"
WORK="/cfs/klemming/projects/supr/naiss2025-23-666/Saida/Scripts/vcf_output/PANGENOME/graph_pangenome"
VCF=$WORK/vcf.noDUPTANDEM_keepINV.vcf.gz

mkdir -p $WORK/03_graph $WORK/04_index $WORK/05_logs $WORK/06_qc

echo "== Sanity checks =="
test -f "$REF" || (echo "Missing REF: $REF" && exit 1)
test -f "$VCF" || (echo "Missing VCF: $VCF" && exit 1)
test -f "${REF}.fai" || (echo "Missing FASTA index: ${REF}.fai" && exit 1)
test -f "${VCF}.tbi" || (echo "Missing VCF index: ${VCF}.tbi" && exit 1)

echo "== Construct graph (SV graph) =="
# -a: keep allele paths (useful if later genotyping known variants)
# -S: allow symbolic SV parsing where possible

# Note: For some symbolic types vg may skip; that's OK for a 'clean' SV sequence graph.
vg construct -r "$REF" -v "$VCF" -a -S > $WORK/03_graph/salix.sv.vg

echo "== Graph stats =="
vg stats -z $WORK/03_graph/salix.sv.vg > $WORK/06_qc/salix.sv.vg.stats.txt || true

echo "== Prune graph for mapping index (GCSA) =="
vg prune -r $WORK/03_graph/salix.sv.vg > $WORK/03_graph/salix.sv.pruned.vg

echo "== Build XG index (unpruned; keep allele paths) =="
vg index -x $WORK/04_index/salix.sv.xg -L $WORK/03_graph/salix.sv.vg

echo "== Build GCSA index (pruned) =="
vg index -g $WORK/04_index/salix.sv.gcsa -k 16 $WORK/03_graph/salix.sv.pruned.vg

echo "== Done =="
ls -lh $WORK/03_graph $WORK/04_index

