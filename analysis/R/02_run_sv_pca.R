#!/usr/bin/env Rscript

# ============================================================
# Structural Variant PCA from presence/absence genotypes
#
# Author: Saida Sharifova
# Project: Structural Variant Analysis and Pangenome Construction
#          in Salix viminalis
# Description:
#   This script reads an SV genotype table, converts genotypes
#   to presence/absence, performs PCA across samples, writes
#   PCA outputs, and produces a labeled PCA plot.
#
# Input files:
#   sv_GT.tsv
#   samples.txt
#   meta.tsv
#
# Output files:
#   SV_PCA_scores.tsv
#   SV_PCA_variance.tsv
#   SV_PCA_plot.pdf
#
# Date: 2026-03-24
# ============================================================

# ------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# ------------------------------------------------------------
# Input / output file names
# ------------------------------------------------------------
genotype_file <- "sv_GT.tsv"
sample_file   <- "samples.txt"
metadata_file <- "meta.tsv"

scores_file   <- "SV_PCA_scores.tsv"
variance_file <- "SV_PCA_variance.tsv"
plot_file     <- "SV_PCA_plot.pdf"

# ------------------------------------------------------------
# Step 1: Read SV genotype table
# Expected format:
#   CHROM  POS  ID  sample1_GT  sample2_GT ...
# ------------------------------------------------------------
sv_table <- read.table(
  genotype_file,
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Create unique SV identifiers
sv_ids <- paste(sv_table[, 1], sv_table[, 2], sv_table[, 3], sep = ":")

# Extract genotype columns only
genotype_table <- sv_table[, -(1:3), drop = FALSE]

# ------------------------------------------------------------
# Step 2: Convert genotypes to presence/absence
#
# Coding:
#   0/0, ./., .|. -> 0
#   0/1, 1/0, 1/1, 0|1, 1|0, 1|1 -> 1
#
# Note:
#   Even if current data are unphased, phased versions are also
#   handled here to make the script more robust.
# ------------------------------------------------------------
genotype_to_pa <- function(g) {
  out <- rep(0L, length(g))
  out[g %in% c("0/1", "1/0", "1/1", "0|1", "1|0", "1|1")] <- 1L
  return(out)
}

# Apply conversion to each sample column
# Current orientation: rows = SVs, columns = samples
sv_by_sample_matrix <- as.matrix(
  as.data.frame(lapply(genotype_table, genotype_to_pa))
)

# Transpose to rows = samples, columns = SVs
sample_by_sv_matrix <- t(sv_by_sample_matrix)

# ------------------------------------------------------------
# Step 3: Assign sample names
# ------------------------------------------------------------
sample_names <- readLines(sample_file)

if (length(sample_names) != nrow(sample_by_sv_matrix)) {
  stop("Number of sample names does not match the number of genotype columns.")
}

rownames(sample_by_sv_matrix) <- sample_names
colnames(sample_by_sv_matrix) <- sv_ids

# ------------------------------------------------------------
# Step 4: Remove invariant SVs
# These have zero variance across samples and should be removed
# before PCA with scaling.
# ------------------------------------------------------------
variable_sv <- apply(sample_by_sv_matrix, 2, var) > 0
sample_by_sv_matrix_filtered <- sample_by_sv_matrix[, variable_sv, drop = FALSE]

cat("Total SVs:", ncol(sample_by_sv_matrix), "\n")
cat("Variable SVs retained for PCA:", ncol(sample_by_sv_matrix_filtered), "\n")
cat("Total samples:", nrow(sample_by_sv_matrix_filtered), "\n")

# ------------------------------------------------------------
# Step 5: Run PCA
# ------------------------------------------------------------
sv_pca <- prcomp(
  sample_by_sv_matrix_filtered,
  center = TRUE,
  scale. = TRUE
)

# ------------------------------------------------------------
# Step 6: Write PCA scores
# ------------------------------------------------------------
pca_scores <- data.frame(
  sample = rownames(sv_pca$x),
  PC1 = sv_pca$x[, 1],
  PC2 = sv_pca$x[, 2],
  stringsAsFactors = FALSE
)

write.table(
  pca_scores,
  file = scores_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Step 7: Write variance explained
# ------------------------------------------------------------
variance_explained <- (sv_pca$sdev^2) / sum(sv_pca$sdev^2)

pca_variance <- data.frame(
  PC = seq_along(variance_explained),
  VarExplained = variance_explained
)

write.table(
  pca_variance,
  file = variance_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# Step 8: Read metadata and merge with PCA scores
# Expected metadata columns:
#   sample   group
# ------------------------------------------------------------
metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

plot_data <- merge(pca_scores, metadata, by = "sample", all.x = TRUE)

# Stop if any samples are missing metadata
if (any(is.na(plot_data$group))) {
  missing_samples <- plot_data$sample[is.na(plot_data$group)]
  stop(
    paste(
      "The following samples are missing group information in meta.tsv:\n",
      paste(missing_samples, collapse = "\n")
    )
  )
}

# ------------------------------------------------------------
# Step 9: Recode group labels
# ------------------------------------------------------------
plot_data$group <- as.character(plot_data$group)

plot_data$group <- recode(
  plot_data$group,
  "S" = "Sweden",
  "U" = "United Kingdom",
  "W" = "Western Europe",
  "E" = "Eastern Europe",
  "R" = "Russia",
  "Breeding_populations" = "Breeding populations",
  .default = plot_data$group
)

plot_data$group <- factor(
  plot_data$group,
  levels = c(
    "Sweden",
    "United Kingdom",
    "Western Europe",
    "Eastern Europe",
    "Russia",
    "Breeding populations"
  )
)

# ------------------------------------------------------------
# Step 10: Create shorter sample labels
# Example:
#   TK-2744-1001_S173 -> 1001
# If pattern does not match, keep original sample name.
# ------------------------------------------------------------
plot_data$label <- ifelse(
  grepl("^TK-[^-]+-[0-9]+_.*$", plot_data$sample),
  sub("^TK-[^-]+-([0-9]+)_.*$", "\\1", plot_data$sample),
  plot_data$sample
)

# ------------------------------------------------------------
# Step 11: Axis labels with explained variance
# ------------------------------------------------------------
pc1_percent <- round(100 * pca_variance$VarExplained[pca_variance$PC == 1], 2)
pc2_percent <- round(100 * pca_variance$VarExplained[pca_variance$PC == 2], 2)

# ------------------------------------------------------------
# Step 12: Keep valid rows for plotting
# ------------------------------------------------------------
plot_data_clean <- plot_data %>%
  filter(is.finite(PC1), is.finite(PC2), !is.na(group))

cat("\nSamples per group:\n")
print(table(plot_data_clean$group))

# ------------------------------------------------------------
# Step 13: PCA plot
# ------------------------------------------------------------
pca_plot <- ggplot(plot_data_clean, aes(x = PC1, y = PC2, color = group)) +
  
  stat_ellipse(
    aes(fill = group),
    type = "t",
    level = 0.68,
    geom = "polygon",
    alpha = 0.08,
    color = NA,
    na.rm = TRUE
  ) +
  
  stat_ellipse(
    type = "t",
    level = 0.68,
    linewidth = 0.35,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  
  geom_point(
    shape = 17,
    size = 2.1,
    alpha = 0.9,
    na.rm = TRUE
  ) +
  
  geom_text_repel(
    aes(label = label),
    size = 1.7,
    box.padding = 0.12,
    point.padding = 0.10,
    max.overlaps = Inf,
    seed = 1,
    min.segment.length = 0,
    segment.size = 0.08,
    segment.alpha = 0.25,
    force = 0.6,
    max.time = 20,
    na.rm = TRUE
  ) +
  
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  coord_equal() +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "PCA of Structural Variants (presence/absence)",
    x = paste0("PC1 (", pc1_percent, "%)"),
    y = paste0("PC2 (", pc2_percent, "%)"),
    color = "Origin",
    fill = "Origin"
  )

# ------------------------------------------------------------
# Step 14: Save plot
# ------------------------------------------------------------
ggsave(
  filename = plot_file,
  plot = pca_plot,
  width = 9,
  height = 7,
  units = "in"
)

cat("\nOutput files written:\n")
cat(" -", scores_file, "\n")
cat(" -", variance_file, "\n")
cat(" -", plot_file, "\n")
