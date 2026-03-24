#!/usr/bin/env Rscript

# ============================================================
# 04_graph_sv_population_analysis.R
#
# Purpose:
#   Population-level analysis of graph-based SVs (already filtered ≥50 bp)
#
# Generates:
#     1. SV counts per sample by population
#     2. Population-specific AF distribution
#     3. Population-specific MAF distribution
#
# Input:
#   graph_sv_table.tsv
#   sample_population.tsv
#
# Input preparation:
#   bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
#     salix_graph_75samples.merged.PASS.SV50.vcf.gz > graph_sv_table.tsv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
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
# Step 1: Read SV table
# ------------------------------------------------------------
sv_raw <- read_tsv(sv_table_file, show_col_types = FALSE)

# Clean column names
colnames(sv_raw) <- colnames(sv_raw) |>
  gsub("^#\\[[0-9]+\\]", "", x = _) |>
  gsub("^\\[[0-9]+\\]", "", x = _) |>
  gsub(":GT$", "", x = _)

# ------------------------------------------------------------
# Step 2: Infer SV type
# ------------------------------------------------------------
sv_annotated <- sv_raw %>%
  mutate(
    SVLEN = nchar(ALT) - nchar(REF),
    SVTYPE = case_when(
      SVLEN > 0 ~ "INS",
      SVLEN < 0 ~ "DEL",
      TRUE ~ "OTHER"
    )
  )

# ------------------------------------------------------------
# Step 3: Identify sample columns
# ------------------------------------------------------------
sample_cols <- colnames(sv_annotated)[grepl("^Sample_", colnames(sv_annotated))]

# ------------------------------------------------------------
# Step 4: Long format
# ------------------------------------------------------------
sv_long <- sv_annotated %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "GT"
  )

# ------------------------------------------------------------
# Step 5: SV counts per sample
# ------------------------------------------------------------
sv_counts <- sv_long %>%
  mutate(
    present = case_when(
      GT %in% c("0/1", "1/0", "1/1") ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  filter(present == 1) %>%
  group_by(sample, SVTYPE) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(sample = str_remove(sample, "^Sample_"))

# ------------------------------------------------------------
# Step 6: Population metadata
# ------------------------------------------------------------
population_df <- read_tsv(population_file, show_col_types = FALSE) %>%
  mutate(sample = str_remove(sample, "^Sample_"))

# ------------------------------------------------------------
# Step 7: Merge
# ------------------------------------------------------------
sv_pop <- sv_counts %>%
  left_join(population_df, by = "sample") %>%
  mutate(population = recode_population(population))

# ------------------------------------------------------------
# Step 8: Labels 
# ------------------------------------------------------------
label_df <- population_df %>%
  mutate(population = recode_population(population)) %>%
  count(population) %>%
  mutate(label = paste0("n = ", n))

# Position labels ABOVE data
y_min <- min(sv_pop$count)

# ------------------------------------------------------------
# Plot 1: SV counts by population
# ------------------------------------------------------------
p <- ggplot(sv_pop, aes(x = population, y = count, fill = SVTYPE)) +
  
  geom_boxplot(
    outlier.shape = NA,
    width = 0.7,
    position = position_dodge(width = 0.8),
    color = "black"
  ) +
  
  geom_jitter(
    aes(color = SVTYPE),
    size = 1.3,
    alpha = 0.5,
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = 0.8
    )
  ) +
  
  #  LABEL POSITION
  geom_text(
    data = label_df,
    aes(x = population, y = y_min + 800, label = label),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) +
  
  scale_fill_manual(
    values = c("DEL" = "#1B9E77", "INS" = "#D95F02"),
    labels = c("DEL" = "Deletion", "INS" = "Insertion")
  ) +
  
  scale_color_manual(
    values = c("DEL" = "#1B9E77", "INS" = "#D95F02"),
    labels = c("DEL" = "Deletion", "INS" = "Insertion")
  ) +
  
  labs(
    x = "Populations",
    y = "Number of structural variants"
  ) +
  
  coord_cartesian(ylim = c(1200, 12000)) +
  
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 20, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5)
  )

ggsave("SV_counts_by_population.pdf", plot = p, width = 9, height = 5)
ggsave("SV_counts_by_population.png", plot = p, width = 9, height = 5, dpi = 600)

# ------------------------------------------------------------
# Step 9: Population-specific AF
# ------------------------------------------------------------
sv_long_pop <- sv_long %>%
  mutate(sample = str_remove(sample, "^Sample_")) %>%
  left_join(population_df, by = "sample") %>%
  mutate(
    population = recode_population(population),
    allele_count = case_when(
      GT == "0/0" ~ 0,
      GT %in% c("0/1", "1/0") ~ 1,
      GT == "1/1" ~ 2,
      TRUE ~ NA_real_
    )
  )

sv_af_pop <- sv_long_pop %>%
  group_by(ID, SVTYPE, population) %>%
  summarise(
    AC = sum(allele_count, na.rm = TRUE),
    AN = sum(!is.na(allele_count)) * 2,
    AF = AC / AN,
    .groups = "drop"
  ) %>%
  filter(!is.na(population), AN > 0)

p2 <- ggplot(sv_af_pop, aes(x = AF, fill = SVTYPE)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ population, scales = "free_y") +
  scale_fill_manual(values = c("DEL" = "#1B9E77", "INS" = "#D95F02")) +
  theme_classic()

ggsave("AF_by_population.pdf", plot = p2, width = 10, height = 6)
ggsave("AF_by_population.png", plot = p2, width = 9, height = 5, dpi = 600)

# ------------------------------------------------------------
# Step 10: Population-specific MAF
# ------------------------------------------------------------
sv_maf_pop <- sv_af_pop %>%
  mutate(MAF = pmin(AF, 1 - AF))

p3 <- ggplot(sv_maf_pop, aes(x = MAF, fill = SVTYPE)) +
  geom_histogram(bins = 30, color = "white") +
  facet_wrap(~ population, scales = "free_y") +
  scale_fill_manual(values = c("DEL" = "#1B9E77", "INS" = "#D95F02")) +
  theme_classic()

ggsave("MAF_by_population.pdf", plot = p3, width = 10, height = 6)
ggsave("MAF_by_population.png", plot = p3, width = 9, height = 5, dpi = 600)

cat("04_graph_sv_population_analysis.R completed successfully.\n")
