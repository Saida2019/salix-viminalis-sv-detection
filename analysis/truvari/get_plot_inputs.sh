#!/bin/bash

set -euo pipefail

module load bioinfo-tools
module load bcftools
module load htslib


WORKDIR="/path/to/project/remapping/TRUVARI"

LINEAR="${WORKDIR}/salix_linear_PASS_SV50.vcf.gz"
GRAPH="${WORKDIR}/salix_graph_PASS_SV50.vcf.gz"
TRUVARI_OUT="${WORKDIR}/truvari_final"

OUTDIR="${WORKDIR}/plot_inputs"
mkdir -p "$OUTDIR"

# SANITY CHECKS
test -f "$LINEAR" || { echo "ERROR: Missing linear VCF"; exit 1; }
test -f "$GRAPH" || { echo "ERROR: Missing graph VCF"; exit 1; }
test -d "$TRUVARI_OUT" || { echo "ERROR: Missing Truvari output"; exit 1; }

# SV TYPE COUNTS
echo "== SV type counts: linear =="
bcftools query -f '%INFO/SVTYPE\n' "$LINEAR" \
| sort | uniq -c \
> "$OUTDIR/linear_sv_type_counts.txt"

echo "== SV type counts: graph =="
bcftools query -f '%ID\n' "$GRAPH" \
| cut -d':' -f1 \
| sed 's/Manta//' \
| sort | uniq -c \
> "$OUTDIR/graph_sv_type_counts.txt"

# SIZE BIN COUNTS
echo "== Size counts: linear =="
bcftools query -f '%INFO/SVLEN\n' "$LINEAR" \
| awk '{
    d=$1; if(d<0) d=-d;
    if(d>=50 && d<100) print "50-99 bp";
    else if(d>=100 && d<1000) print "100-999 bp";
    else if(d>=1000 && d<10000) print "1-9.9 kb";
    else if(d>=10000) print ">=10 kb";
}' \
| sort | uniq -c \
> "$OUTDIR/linear_size_counts.txt"

echo "== Size counts: graph =="
bcftools query -f '%REF\t%ALT\n' "$GRAPH" \
| awk '{
    d=length($1)-length($2); if(d<0) d=-d;
    if(d>=50 && d<100) print "50-99 bp";
    else if(d>=100 && d<1000) print "100-999 bp";
    else if(d>=1000 && d<10000) print "1-9.9 kb";
    else if(d>=10000) print ">=10 kb";
}' \
| sort | uniq -c \
> "$OUTDIR/graph_size_counts.txt"

# TRUVARI COUNTS
echo "== Truvari counts =="

bcftools view -H "$TRUVARI_OUT/fn.vcf" | wc -l > "$OUTDIR/linear_only.txt"
bcftools view -H "$TRUVARI_OUT/fp.vcf" | wc -l > "$OUTDIR/graph_only.txt"
bcftools view -H "$TRUVARI_OUT/tp-base.vcf" | wc -l > "$OUTDIR/shared.txt"

# TRUVARI METRICS
cp "$TRUVARI_OUT/summary.txt" "$OUTDIR/truvari_summary.txt"

echo "== Done =="
ls -lh "$OUTDIR"
