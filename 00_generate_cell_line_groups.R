#!/usr/bin/env Rscript

################################################################################
# 00. Generate DepMap loss-of-function cell-line groups for POG BLOF genes
#
# Purpose
#   Generate the per-gene DepMap mutant/control group files used by the
#   downstream metaGRETTA RNAi and CRISPR knockout screening workflows.
#
# Starting point
#   Supplemental Table S1B: all biallelic loss-of-function (BLOF) alterations
#   identified in the POG cohort.
#
# For each unique POG BLOF gene, GRETTA::select_cell_lines() is used to assign
# DepMap models to the standard GRETTA alteration groups.
#
# BLOF group hierarchy used in this study
#   GRETTA returns Control, HetDel, T-HetDel, and HomDel groups. For downstream
#   BLOF screening, mutant groups are prioritised as follows:
#
#     1. HomDel
#     2. T-HetDel, used only when HomDel is not sufficiently represented
#     3. HetDel is never used as a BLOF mutant group
#
#   This script retains the original GRETTA group assignments. The HomDel ->
#   T-HetDel fallback is applied downstream after restricting cell lines by
#   assay, disease context, and microsatellite status, because group
#   availability can differ between screening contexts.
#
# Prerequisites
#   1. Supplemental Table S1B with columns:
#        POG_patient_id, POG_sample_id, POG_sample_disease, BLOF_gene,
#        amino_acid_p, mut_data
#   2. GRETTA-formatted DepMap 23Q4 data.
#
# Outputs
#   For each POG BLOF gene:
#     <output_dir>/individual_genes/<gene>/cell_line_groups.csv
#
#   A summary of the GRETTA groups available for each POG BLOF gene:
#     <output_dir>/DepMap_POG_BLOF_cell_line_group_summary.csv
#
# Notes
#   - Group assignment is performed by GRETTA::select_cell_lines().
#   - No minimum cell-line count is imposed in this script.
#   - Screenability is evaluated downstream after applying assay, disease,
#     and MSS/MSI restrictions.
#   - HomDel and T-HetDel are not combined to meet the minimum group size.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(GRETTA)
})

# ------------------------------------------------------------------------------
# User-defined paths
# ------------------------------------------------------------------------------

# Supplemental Table S1B defines the patient-derived BLOF genes to screen.
supp_table_s1b_file <- "path/to/Supplemental_Table_S1B.csv"

# GRETTA-formatted DepMap 23Q4 data.
gretta_data_dir <- "path/to/GRETTA_data/23Q4/data"

# Parent directory used by the downstream screening scripts.
output_dir <- "path/to/metaGRETTA_screen_output"

# ------------------------------------------------------------------------------
# Validate input files
# ------------------------------------------------------------------------------

if (!file.exists(supp_table_s1b_file)) {
  stop(
    "Supplemental Table S1B was not found: ",
    supp_table_s1b_file
  )
}

if (!dir.exists(gretta_data_dir)) {
  stop(
    "GRETTA data directory was not found: ",
    gretta_data_dir
  )
}

# ------------------------------------------------------------------------------
# Load Supplemental Table S1B
# ------------------------------------------------------------------------------

supp_table_s1b <- read_csv(
  supp_table_s1b_file,
  show_col_types = FALSE
)

required_s1b_columns <- c(
  "POG_patient_id",
  "POG_sample_id",
  "POG_sample_disease",
  "BLOF_gene",
  "amino_acid_p",
  "mut_data"
)

missing_s1b_columns <- setdiff(
  required_s1b_columns,
  colnames(supp_table_s1b)
)

if (length(missing_s1b_columns) > 0) {
  stop(
    "Supplemental Table S1B is missing required columns: ",
    paste(missing_s1b_columns, collapse = ", ")
  )
}

# S1B defines the patient-derived genes for which DepMap LOF groups are needed.
pog_blof_genes <- supp_table_s1b %>%
  filter(
    !is.na(BLOF_gene),
    BLOF_gene != ""
  ) %>%
  pull(BLOF_gene) %>%
  unique() %>%
  sort()

message(
  length(pog_blof_genes),
  " unique POG BLOF genes found in Supplemental Table S1B."
)

# ------------------------------------------------------------------------------
# Generate standard GRETTA cell-line groups
# ------------------------------------------------------------------------------

individual_gene_dir <- file.path(
  output_dir,
  "individual_genes"
)

dir.create(
  individual_gene_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

group_summary <- map_dfr(
  seq_along(pog_blof_genes),
  function(i) {

    gene <- pog_blof_genes[[i]]

    gene_dir <- file.path(
      individual_gene_dir,
      gene
    )

    group_file <- file.path(
      gene_dir,
      "cell_line_groups.csv"
    )

    dir.create(
      gene_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    message(
      "[", i, "/", length(pog_blof_genes), "] ",
      gene
    )

    # --------------------------------------------------------------------------
    # Reuse an existing grouping file when available
    # --------------------------------------------------------------------------

    if (file.exists(group_file)) {

      cell_line_groups <- read_csv(
        group_file,
        show_col_types = FALSE
      )

      status <- "existing"

    } else {

      # ------------------------------------------------------------------------
      # Generate standard GRETTA alteration groups
      # ------------------------------------------------------------------------

      cell_line_groups <- tryCatch(
        GRETTA::select_cell_lines(
          input_gene = gene,
          data_dir = gretta_data_dir
        ),
        error = function(e) {

          warning(
            "Could not generate cell-line groups for ",
            gene,
            ": ",
            conditionMessage(e)
          )

          NULL
        }
      )

      if (is.null(cell_line_groups) || nrow(cell_line_groups) == 0) {

        return(
          tibble(
            GeneName = gene,
            Group = NA_character_,
            n_cell_lines = 0L,
            BLOF_group_eligible = FALSE,
            BLOF_priority = NA_integer_,
            status = "no_groups_generated"
          )
        )
      }

      write_csv(
        cell_line_groups,
        group_file
      )

      status <- "generated"
    }

    # --------------------------------------------------------------------------
    # Validate GRETTA output
    # --------------------------------------------------------------------------

    required_group_columns <- c(
      "DepMap_ID",
      "disease",
      "Group"
    )

    missing_group_columns <- setdiff(
      required_group_columns,
      colnames(cell_line_groups)
    )

    if (length(missing_group_columns) > 0) {
      stop(
        "Cell-line grouping for ",
        gene,
        " is missing required columns: ",
        paste(missing_group_columns, collapse = ", ")
      )
    }

    # --------------------------------------------------------------------------
    # Summarise groups and document the BLOF hierarchy
    #
    # IMPORTANT:
    # The hierarchy is documented here but is not resolved here. Whether
    # HomDel has enough cell lines must be determined within the final
    # assay/disease/MS-status context.
    #
    # HomDel   -> preferred BLOF group
    # T-HetDel -> fallback if HomDel has insufficient representation
    # HetDel   -> never used for BLOF screening
    # Control  -> comparison group
    # --------------------------------------------------------------------------

    cell_line_groups %>%
      count(
        Group,
        name = "n_cell_lines"
      ) %>%
      mutate(
        GeneName = gene,

        BLOF_group_eligible = Group %in% c(
          "HomDel",
          "T-HetDel"
        ),

        BLOF_priority = case_when(
          Group == "HomDel" ~ 1L,
          Group == "T-HetDel" ~ 2L,
          TRUE ~ NA_integer_
        ),

        status = status
      ) %>%
      select(
        GeneName,
        Group,
        n_cell_lines,
        BLOF_group_eligible,
        BLOF_priority,
        status
      )
  }
)

# ------------------------------------------------------------------------------
# Save summary
# ------------------------------------------------------------------------------

summary_file <- file.path(
  output_dir,
  "DepMap_POG_BLOF_cell_line_group_summary.csv"
)

write_csv(
  group_summary,
  summary_file
)

# ------------------------------------------------------------------------------
# Report
# ------------------------------------------------------------------------------

message(
  "\nCell-line grouping complete."
)

message(
  "Per-gene grouping files: ",
  individual_gene_dir
)

message(
  "Group summary: ",
  summary_file
)

message(
  "\nDownstream BLOF hierarchy: HomDel -> T-HetDel fallback; HetDel excluded."
)
