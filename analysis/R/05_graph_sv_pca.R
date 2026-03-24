#!/usr/bin/env Rscript

# ============================================================
# 05_graph_sv_pca.R
#
#   Run PCA on graph-based SV genotypes from a genotype table
#
# Input:
#   graph_sv_table.tsv
#   sample_population.tsv
#
# Output:
#   PCA_graph_75_samples.pdf/png
#
# Input preparation:
#   bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
#     salix_graph_75samples.merged.PASS.SV50.vcf.gz > graph_sv_table.tsv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
sv_table_file <- "graph_sv_table.tsv"
population_file <- "sample_population.tsv"

# ------------------------------------------------------------
# Helper: population labels
# ------------------------------------------------------------
recode_population <- function(x) {
  factor(
    x,
    levels = c("SW", "UK", "R", "E", "WE", "Br"),
    labels = c(
      "Sweden",
      "UK",
      "Russia",
      "Eastern Europe",
      "Western Europe",
      "Breeding populations"
    )
  )
}

# ------------------------------------------------------------
# Step 1: Read graph-based SV table
# ------------------------------------------------------------
sv_raw <- read_tsv(sv_table_file, show_col_types = FALSE)

# Clean bcftools-generated column names
colnames(sv_raw) <- colnames(sv_raw) |>
  gsub("^#\\[[0-9]+\\]", "", x = _) |>
  gsub("^\\[[0-9]+\\]", "", x = _) |>
  gsub(":GT$", "", x = _)

# ------------------------------------------------------------
# Step 2: Add inferred SV annotations
# ------------------------------------------------------------
sv_annotated <- sv_raw %>%
  mutate(
    ref_len = nchar(REF),
    alt_len = nchar(ALT),
    SVLEN = alt_len - ref_len,
    SVTYPE = case_when(
      SVLEN > 0 ~ "INS",
      SVLEN < 0 ~ "DEL",
      TRUE ~ "OTHER"
    ),
    var_id = paste(CHROM, POS, REF, ALT, sep = "_")
  )

cat("Rows before deduplication:", nrow(sv_annotated), "\n")
cat("Distinct var_id before deduplication:", n_distinct(sv_annotated$var_id), "\n")
cat("Duplicated var_id before deduplication:", sum(duplicated(sv_annotated$var_id)), "\n")

# ------------------------------------------------------------
# Step 3: Remove duplicated variants
# ------------------------------------------------------------
sv_annotated <- sv_annotated %>%
  distinct(var_id, .keep_all = TRUE)

cat("Rows after deduplication:", nrow(sv_annotated), "\n")
cat("Distinct var_id after deduplication:", n_distinct(sv_annotated$var_id), "\n")

# ------------------------------------------------------------
# Step 4: Identify sample genotype columns
# ------------------------------------------------------------
sample_cols <- colnames(sv_annotated)[grepl("^Sample_", colnames(sv_annotated))]

cat("Number of sample columns:", length(sample_cols), "\n")
print(head(sample_cols))

# ------------------------------------------------------------
# Step 5: Convert to long format
# ------------------------------------------------------------
sv_long_pca <- sv_annotated %>%
  select(var_id, all_of(sample_cols)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "GT"
  )

print(dim(sv_long_pca))
print(table(sv_long_pca$GT, useNA = "ifany"))

# ------------------------------------------------------------
# Step 6: Convert genotypes to dosage
#   0/0 = 0
#   0/1 or 1/0 = 1
#   1/1 = 2
#   other or missing = NA
# ------------------------------------------------------------
sv_long_pca <- sv_long_pca %>%
  mutate(
    dosage = case_when(
      GT == "0/0" ~ 0,
      GT %in% c("0/1", "1/0") ~ 1,
      GT == "1/1" ~ 2,
      TRUE ~ NA_real_
    )
  )

print(table(sv_long_pca$dosage, useNA = "ifany"))

# ------------------------------------------------------------
# Step 7: Check duplicate sample-variant combinations
# ------------------------------------------------------------
dup_check <- sv_long_pca %>%
  summarise(n = n(), .by = c(sample, var_id)) %>%
  filter(n > 1)

cat("Duplicated sample-variant combinations:", nrow(dup_check), "\n")

# ------------------------------------------------------------
# Step 8: Build wide genotype matrix
# ------------------------------------------------------------
sv_pca_wide <- sv_long_pca %>%
  select(sample, var_id, dosage) %>%
  pivot_wider(
    names_from = var_id,
    values_from = dosage
  )

print(dim(sv_pca_wide))

sample_ids <- sv_pca_wide$sample

geno_mat <- sv_pca_wide %>%
  select(-sample) %>%
  as.data.frame()

cat("Missing values before imputation:", sum(is.na(geno_mat)), "\n")

# ------------------------------------------------------------
# Step 9: Mean-impute missing values by variant
# ------------------------------------------------------------
geno_mat <- geno_mat %>%
  mutate(across(
    everything(),
    ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)
  ))

cat("Missing values after imputation:", sum(is.na(geno_mat)), "\n")

# ------------------------------------------------------------
# Step 10: Remove invariant variants
# ------------------------------------------------------------
variant_variance <- apply(geno_mat, 2, var)
cat("Invariant variants:", sum(variant_variance == 0), "\n")

geno_mat <- geno_mat[, variant_variance > 0, drop = FALSE]

cat("Final PCA matrix dimensions:\n")
print(dim(geno_mat))

# ------------------------------------------------------------
# Step 11: Run PCA
# ------------------------------------------------------------
geno_mat <- as.matrix(geno_mat)
rownames(geno_mat) <- sample_ids

pca_res <- prcomp(geno_mat, center = TRUE, scale. = TRUE)

pve <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
cat("Variance explained by first 10 PCs:\n")
print(round(pve[1:10], 2))

# ------------------------------------------------------------
# Step 12: Read population metadata
# ------------------------------------------------------------
population_df <- read_tsv(population_file, show_col_types = FALSE) %>%
  mutate(sample = str_remove(sample, "^Sample_"))

# ------------------------------------------------------------
# Step 13: Build PCA plotting table
# ------------------------------------------------------------
pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  mutate(sample = str_remove(sample, "^Sample_")) %>%
  left_join(population_df, by = "sample") %>%
  mutate(
    population = recode_population(population),
    sample_label = sample,
    PC1 = -PC1
  )

print(table(pca_df$population, useNA = "ifany"))

# ------------------------------------------------------------
# Step 14: Select samples to label
# ------------------------------------------------------------
label_samples <- pca_df %>%
  filter(
    population %in% c("Russia", "Breeding populations") |
      sample_label %in% c("TK-2744-73", "TK-2744-310", "TK-2744-308")
  )

# ------------------------------------------------------------
# Step 15: PCA plot
# ------------------------------------------------------------
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = population, shape = population)) +
  stat_ellipse(
    type = "norm",
    level = 0.68,
    linewidth = 0.8,
    alpha = 0.10,
    show.legend = FALSE
  ) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(
    data = label_samples,
    aes(label = sample_label),
    size = 2.5,
    max.overlaps = 100,
    box.padding = 0.35,
    point.padding = 0.25,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Sweden" = "#1F78B4",
    "UK" = "#E66100",
    "Russia" = "#6A3D9A",
    "Eastern Europe" = "#D81B60",
    "Western Europe" = "#33A02C",
    "Breeding populations" = "#C49A00"
  )) +
  scale_shape_manual(values = c(
    "Sweden" = 16,
    "UK" = 17,
    "Russia" = 15,
    "Eastern Europe" = 18,
    "Western Europe" = 16,
    "Breeding populations" = 17
  )) +
  labs(
    x = paste0("PC1 (", round(pve[1], 1), "%)"),
    y = paste0("PC2 (", round(pve[2], 1), "%)")
  ) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5)
  )

ggsave("PCA_graph_75_samples.pdf", plot = p_pca, width = 9, height = 5)
ggsave("PCA_graph_75_samples.png", plot = p_pca, width = 9, height = 5, dpi = 600)

cat("05_graph_sv_pca.R completed successfully.\n")
