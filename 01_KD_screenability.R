#!/usr/bin/env Rscript

################################################################################
# 01. Identify screenable POG BLOF contexts for RNAi dependency screening
#
# Starting point
#   Supplemental Table S1B: all biallelic loss-of-function (BLOF) alterations
#   identified in the POG cohort.
#
# Purpose
#   For each POG BLOF gene represented in DepMap LOF groups, determine which
#   gene x LOF-group x disease x microsatellite-status contexts contain at
#   least the minimum number of mutant and control cell lines with RNAi data.
#
# Required columns in Supplemental Table S1B:
#   POG_patient_id, POG_sample_id, POG_sample_disease, BLOF_gene,
#   amino_acid_p, mut_data
#
# Required per-gene file:
#   <lof_context_dir>/<gene>/cell_line_groups.csv
# with columns:
#   DepMap_ID, disease, Group
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(GRETTA)
})

# ------------------------------------------------------------------------------
# User-defined paths and parameters
# ------------------------------------------------------------------------------

supp_table_s1b_file <- "path/to/Supplemental_Table_S1B.csv"

# GRETTA-formatted DepMap 23Q4 data.
gretta_data_dir <- "path/to/GRETTA_data/23Q4/data"

# Directory containing per-gene `cell_line_groups.csv` files.
lof_context_dir <- "path/to/individual_genes"

# Output directory for the RNAi workflow.
output_dir <- "results/rnai"

min_group_size <- 3L

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------

supp_table_s1b <- read_csv(supp_table_s1b_file, show_col_types = FALSE)

required_s1b_columns <- c(
  "POG_patient_id",
  "POG_sample_id",
  "POG_sample_disease",
  "BLOF_gene",
  "amino_acid_p",
  "mut_data"
)

missing_s1b_columns <- setdiff(required_s1b_columns, colnames(supp_table_s1b))
if (length(missing_s1b_columns) > 0) {
  stop(
    "Supplemental Table S1B is missing required columns: ",
    paste(missing_s1b_columns, collapse = ", ")
  )
}

load(file.path(gretta_data_dir, "rnai_long.rda"))

ms_status_df <- read_csv(
  file.path(gretta_data_dir, "inferred_MS_status.csv"),
  show_col_types = FALSE
) %>%
  filter(Inferred_MS_status %in% c("inferred_MSS", "inferred_MSI")) %>%
  distinct(DepMap_ID, Inferred_MS_status)

stopifnot("DepMap_ID" %in% colnames(rnai_long))

# POG BLOF genes define the query-gene universe for the publication analysis.
pog_blof_genes <- supp_table_s1b %>%
  filter(!is.na(BLOF_gene), BLOF_gene != "") %>%
  pull(BLOF_gene) %>%
  unique() %>%
  sort()

message(
  length(pog_blof_genes),
  " unique POG BLOF genes found in Supplemental Table S1B."
)

# Disease contexts available in GRETTA.
disease_list <- GRETTA::list_cancer_types(data_dir = gretta_data_dir)
disease_list <- disease_list[
  !is.na(disease_list) &
    !disease_list %in% c("Normal", "Embryonal", "Other")
]
disease_list <- c("pancan", disease_list)

# Only models with RNAi data and a known MSS/MSI assignment are eligible.
eligible_rnai_ids <- intersect(
  unique(rnai_long$DepMap_ID),
  unique(ms_status_df$DepMap_ID)
)

# ------------------------------------------------------------------------------
# Identify screenable contexts
# ------------------------------------------------------------------------------

screenable_contexts <- map_dfr(pog_blof_genes, function(gene) {

  group_file <- file.path(
    lof_context_dir,
    gene,
    "cell_line_groups.csv"
  )

  if (!file.exists(group_file)) {
    return(tibble())
  }

  cell_lines <- read_csv(group_file, show_col_types = FALSE)

  required_group_columns <- c("DepMap_ID", "disease", "Group")
  missing_group_columns <- setdiff(required_group_columns, colnames(cell_lines))

  if (length(missing_group_columns) > 0) {
    stop(
      "Missing required columns in ", group_file, ": ",
      paste(missing_group_columns, collapse = ", ")
    )
  }

  # Apply the BLOF hierarchy independently within each disease x MS-status
  # context:
  #   1. HomDel is preferred when >= min_group_size cell lines are available.
  #   2. T-HetDel is used only when HomDel is insufficient and T-HetDel has
  #      >= min_group_size cell lines.
  #   3. HetDel is never used.
  # HomDel and T-HetDel are not pooled to satisfy the minimum group size.
  crossing(
    disease = disease_list,
    MS_status = c("MSS", "MSI")
  ) %>%
    mutate(
      counts = map2(
        disease,
        MS_status,
        function(disease_name, ms_label) {

          inferred_status <- if_else(
            ms_label == "MSS",
            "inferred_MSS",
            "inferred_MSI"
          )

          ms_ids <- ms_status_df %>%
            filter(Inferred_MS_status == inferred_status) %>%
            pull(DepMap_ID)

          context_lines <- cell_lines %>%
            filter(
              Group %in% c("Control", "HomDel", "T-HetDel"),
              DepMap_ID %in% eligible_rnai_ids,
              DepMap_ID %in% ms_ids
            )

          if (disease_name != "pancan") {
            context_lines <- context_lines %>%
              filter(disease == disease_name)
          }

          n_control <- context_lines %>%
            filter(Group == "Control") %>%
            pull(DepMap_ID) %>%
            unique() %>%
            length()

          n_homdel <- context_lines %>%
            filter(Group == "HomDel") %>%
            pull(DepMap_ID) %>%
            unique() %>%
            length()

          n_thet_del <- context_lines %>%
            filter(Group == "T-HetDel") %>%
            pull(DepMap_ID) %>%
            unique() %>%
            length()

          mut_group <- case_when(
            n_homdel >= min_group_size ~ "HomDel",
            n_thet_del >= min_group_size ~ "T-HetDel",
            TRUE ~ NA_character_
          )

          n_mutant <- case_when(
            mut_group == "HomDel" ~ n_homdel,
            mut_group == "T-HetDel" ~ n_thet_del,
            TRUE ~ 0L
          )

          tibble(
            mut_group = mut_group,
            n_control = n_control,
            n_homdel = n_homdel,
            n_thet_del = n_thet_del,
            n_mutant = n_mutant
          )
        }
      )
    ) %>%
    unnest(counts) %>%
    mutate(
      GeneName = gene,
      is_screenable = (
        n_control >= min_group_size &
          !is.na(mut_group) &
          n_mutant >= min_group_size
      )
    ) %>%
    select(
      GeneName,
      mut_group,
      disease,
      MS_status,
      n_control,
      n_homdel,
      n_thet_del,
      n_mutant,
      is_screenable
    )
})

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

screenable_file <- file.path(
  output_dir,
  "RNAi_screenable_POG_BLOF_contexts.csv"
)

write_csv(screenable_contexts, screenable_file)

message(
  sum(screenable_contexts$is_screenable),
  " screenable RNAi contexts identified."
)
message("Results written to: ", screenable_file)
