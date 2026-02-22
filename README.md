# Structural Variant Detection and Graph Pangenome Construction in *Salix viminalis*

Bioinformatics workflows and scripts for preprocessing, structural variant (SV) detection, and VG-based pangenome analysis in *Salix viminalis*.

---

## Repository Structure

### scripts/preprocessing/
Standalone shell scripts for raw read preprocessing:

- FastQC quality control
- MultiQC summary reports
- BWA mapping
- Picard read group assignment
- Picard duplicate marking

These scripts are designed to run on HPC systems using SLURM and environment modules.

---

### workflows/sv_nextflow/
Nextflow pipeline for structural variant detection.

Contains:
- `main.nf` — workflow definition
- `nextflow.config` — configuration file
- `modules/` — Nextflow workflow modules

Designed for execution on SLURM clusters.

---

### workflows/vg_pangenome/
Shell-based workflow for VG pangenome construction and analysis.

Contains stepwise scripts for:
- Graph construction
- Graph indexing
- Read mapping
- Variant calling

---

## Execution Environment

All workflows are designed for HPC clusters using:

- SLURM scheduler
- Environment modules (e.g. `module load vg`, `bwa`, `picard`, etc.)

---

## Author

Saida Sharifova
