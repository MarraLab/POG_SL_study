#!/usr/bin/env Rscript

################################################################################
# CRISPR knockout synthetic-lethality screening by microsatellite status
#
# Purpose
#   Identify screenable loss-of-function (LOF) contexts and run GRETTA CRISPR
#   knockout (KO) screens separately in microsatellite-stable (MSS) and
#   microsatellite-instability (MSI) DepMap cell lines.
#
# This script consolidates the original MSS and MSI workflows into one
# publication-facing script. The same analysis is run for each microsatellite
# group by changing only the cell-line subset.
#
# Starting point
#   Supplemental Table S1B: all biallelic loss-of-function (BLOF) alterations
#   identified in the POG cohort.
#
# Prerequisites
#   1. Supplemental Table S1B with columns:
#        POG_patient_id, POG_sample_id, POG_sample_disease, BLOF_gene,
#        amino_acid_p, mut_data
#   2. GRETTA-formatted DepMap 23Q4 data.
#   3. A CSV assigning each DepMap model to inferred MSS/MSI status.
#   4. Per-gene `cell_line_groups.csv` files generated during DepMap
#      LOF-context definition. Each file must contain:
#        DepMap_ID, disease, Group
#      where "Control" denotes the comparison group and non-control groups
#      denote LOF contexts.
#   5. Open Targets tractability annotations (v23.02 used in the study).
#
# Outputs
#   For each MSS/MSI subset:
#     - a table of screenable gene / LOF-group / disease contexts
#     - individual GRETTA screen results
#     - Tier 1/2/3 combined results
#     - optional screen plots
#
# Notes
#   - A screen requires >= 3 mutant and >= 3 control cell lines.
#   - "pancan" includes all eligible cancer cell lines in the relevant
#     microsatellite-status subset.
#   - HomDel is preferred as the BLOF mutant group; T-HetDel is used only when
#     HomDel is insufficient. HetDel is never used. HomDel and T-HetDel are
#     not pooled to meet the minimum group size.
#   - The tier definitions below reproduce the CRISPR KO criteria used in the
#     study.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(GRETTA)
  library(janitor)
})

# ------------------------------------------------------------------------------
# User-defined paths and analysis parameters
# ------------------------------------------------------------------------------

# Supplemental Table S1B defines the POG BLOF genes to screen.
supp_table_s1b_file <- "path/to/Supplemental_Table_S1B.csv"

# Directory containing GRETTA-formatted DepMap data.
gretta_data_dir <- "path/to/GRETTA_data/23Q4/data"

# Directory containing the per-gene `cell_line_groups.csv` files generated
# by Script 00. Point this to:
#   <script_00_output_dir>/individual_genes
lof_context_dir <- "path/to/metaGRETTA_screen_output/individual_genes"

# Output directory for the KO screening workflow.
output_dir <- "results/ko"

# Microsatellite-status file with columns:
#   DepMap_ID, Inferred_MS_status
ms_status_file <- file.path(gretta_data_dir, "inferred_MS_status.csv")

# Open Targets tractability table used for target annotation.
tractability_file <- "path/to/tractability_v23-02.tsv"

# Microsatellite subsets to analyse.
ms_statuses <- c(
  MSS = "inferred_MSS",
  MSI = "inferred_MSI"
)

# Analysis parameters.
min_group_size <- 3
n_perm <- 10000
core_num <- 12
make_plots <- TRUE

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

required_files <- c(
  supp_table_s1b_file,
  file.path(gretta_data_dir, "dep.rda"),
  ms_status_file,
  tractability_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input files were not found:\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load Supplemental Table S1B and DepMap data
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

# Supplemental Table S1B defines the patient-derived BLOF genes considered for
# functional screening.
pog_blof_genes <- supp_table_s1b %>%
  filter(!is.na(BLOF_gene), BLOF_gene != "") %>%
  pull(BLOF_gene) %>%
  unique() %>%
  sort()

message(
  length(pog_blof_genes),
  " unique POG BLOF genes found in Supplemental Table S1B."
)

# CRISPR dependency data are used to restrict screening to DepMap models with
# KO dependency measurements.
load(file.path(gretta_data_dir, "dep.rda"))

ms_status_df <- read_csv(
  ms_status_file,
  show_col_types = FALSE
)

required_ms_columns <- c("DepMap_ID", "Inferred_MS_status")
if (!all(required_ms_columns %in% colnames(ms_status_df))) {
  stop(
    "`ms_status_file` must contain: ",
    paste(required_ms_columns, collapse = ", ")
  )
}

ms_status_df <- ms_status_df %>%
  filter(Inferred_MS_status %in% c("inferred_MSS", "inferred_MSI")) %>%
  distinct(DepMap_ID, Inferred_MS_status)

# Cancer types available in GRETTA.
disease_list <- GRETTA::list_cancer_types(data_dir = gretta_data_dir)
disease_list <- disease_list[
  !is.na(disease_list) &
    !disease_list %in% c("Normal", "Embryonal", "Other")
]
disease_list <- c("pancan", disease_list)

# ------------------------------------------------------------------------------
# Open Targets tractability annotation
# ------------------------------------------------------------------------------

opentarget_df <- read_tsv(
  tractability_file,
  show_col_types = FALSE
) %>%
  janitor::clean_names()

compound_tractability <- opentarget_df %>%
  select(
    symbol,
    protein_names,
    top_bucket_sm,
    category_sm,
    top_bucket_ab,
    category_ab,
    top_bucket_protac,
    category_protac,
    top_bucket_othercl,
    category_othercl
  ) %>%
  mutate(
    Drug_group = case_when(
      # WRN was manually assigned to Group 1 in the original analysis.
      symbol == "WRN" ~ "Group_1",

      category_sm == "Clinical_Precedence_sm" ~ "Group_1",
      category_ab == "Clinical_Precedence_ab" ~ "Group_1",
      category_protac == "Clinical_Precedence_protac" ~ "Group_1",
      category_othercl == "Clinical_Precedence_othercl" ~ "Group_1",

      category_sm == "Discovery_Precedence_sm" ~ "Group_2",
      category_ab == "Predicted_Tractable_ab_High_confidence_ab" ~ "Group_2",
      category_protac == "Literature_Precedence_protac" ~ "Group_2",

      category_sm == "Predicted_Tractable_sm" ~ "Group_3",
      category_ab == "Predicted_Tractable_ab_Medium_to_low_confidence_ab" ~ "Group_3",
      category_protac == "Discovery_Opportunity_protac" ~ "Group_3",

      TRUE ~ NA_character_
    )
  )

group_1_genes <- compound_tractability %>%
  filter(Drug_group == "Group_1") %>%
  pull(symbol)

group_2_genes <- compound_tractability %>%
  filter(Drug_group == "Group_2") %>%
  pull(symbol)

group_3_genes <- compound_tractability %>%
  filter(Drug_group == "Group_3") %>%
  pull(symbol)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

safe_name <- function(x) {
  str_replace_all(x, "\\s|/", "_")
}


get_screen_ids <- function(cell_lines, disease, mutant_group, eligible_ids) {

  context_lines <- cell_lines %>%
    filter(
      Group %in% c("Control", mutant_group),
      DepMap_ID %in% eligible_ids
    )

  if (disease != "pancan") {
    context_lines <- context_lines %>%
      filter(.data$disease == disease)
  }

  list(
    mutant_id = context_lines %>%
      filter(Group == mutant_group) %>%
      pull(DepMap_ID) %>%
      unique(),

    control_id = context_lines %>%
      filter(Group == "Control") %>%
      pull(DepMap_ID) %>%
      unique()
  )
}


identify_screenable_contexts <- function(
    lof_genes,
    eligible_ids,
    disease_list,
    lof_context_dir,
    min_group_size = 3) {

  message("Identifying screenable LOF contexts...")

  screenable_contexts <- map_dfr(lof_genes, function(gene) {

    group_file <- file.path(
      lof_context_dir,
      gene,
      "cell_line_groups.csv"
    )

    if (!file.exists(group_file)) {
      return(tibble())
    }

    cell_lines <- read_csv(group_file, show_col_types = FALSE)

    required_columns <- c("DepMap_ID", "disease", "Group")
    if (!all(required_columns %in% colnames(cell_lines))) {
      stop(
        "Missing required columns in ", group_file, ": ",
        paste(setdiff(required_columns, colnames(cell_lines)), collapse = ", ")
      )
    }

    # Apply the BLOF hierarchy independently within each disease context:
    #   1. HomDel is preferred when >= min_group_size cell lines are available.
    #   2. T-HetDel is used only when HomDel is insufficient and T-HetDel has
    #      >= min_group_size cell lines.
    #   3. HetDel is never used.
    # HomDel and T-HetDel are not pooled to satisfy the minimum group size.
    map_dfr(disease_list, function(d) {

      context_lines <- cell_lines %>%
        filter(
          Group %in% c("Control", "HomDel", "T-HetDel"),
          DepMap_ID %in% eligible_ids
        )

      if (d != "pancan") {
        context_lines <- context_lines %>%
          filter(.data$disease == d)
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
        GeneName = gene,
        mut_group = mut_group,
        disease = d,
        n_control = n_control,
        n_homdel = n_homdel,
        n_thet_del = n_thet_del,
        n_mutant = n_mutant,
        is_screenable = (
          n_control >= min_group_size &
            !is.na(mut_group) &
            n_mutant >= min_group_size
        )
      )
    })
  })

  screenable_contexts
}


annotate_ko_results <- function(screen_results_permed) {

  screen_results_permed %>%
    mutate(
      GI_type = if_else(
        Mutant_median > Control_median,
        "SL",
        "AL"
      ),

      # CRISPR KO tier definitions used in the study.
      GI_tier = case_when(
        Pval_perm_adj < 0.01 &
          abs(log2FC_by_median) > 2 &
          (Control_median > 0.5 | Mutant_median > 0.5) ~ "Tier_1",

        Pval_perm_adj < 0.05 &
          abs(log2FC_by_median) > 2 &
          (Control_median > 0.5 | Mutant_median > 0.5) ~ "Tier_2",

        Pval_perm_adj < 0.05 &
          abs(log2FC_by_median) > 2 ~ "Tier_3",

        TRUE ~ NA_character_
      ),

      # This reproduces the original analysis, in which genes without an
      # Open Targets Group 1/2 annotation were assigned to Group 3.
      Drug_group = case_when(
        GeneNames %in% group_1_genes ~ "Group_1",
        GeneNames %in% group_2_genes ~ "Group_2",
        GeneNames %in% group_3_genes ~ "Group_3",
        TRUE ~ "Group_3"
      ),

      GI_tier = factor(
        GI_tier,
        levels = c("Tier_1", "Tier_2", "Tier_3")
      ),

      Drug_group = factor(
        Drug_group,
        levels = c("Group_1", "Group_2", "Group_3")
      ),

      Priority_class = case_when(
        GI_tier == "Tier_1" & Drug_group == "Group_1" ~ "Class_A",
        GI_tier == "Tier_2" & Drug_group == "Group_1" ~ "Class_B",
        GI_tier == "Tier_1" ~ "Class_C",
        GI_tier %in% c("Tier_1", "Tier_2") ~ "Class_D",
        TRUE ~ NA_character_
      ),

      Priority_class = factor(
        Priority_class,
        levels = c("Class_A", "Class_B", "Class_C", "Class_D")
      )
    )
}


run_one_ko_screen <- function(
    gene,
    mutant_group,
    disease,
    eligible_ids,
    ms_label,
    lof_context_dir,
    output_dir,
    gretta_data_dir,
    min_group_size,
    n_perm,
    core_num,
    make_plots = TRUE) {

  group_file <- file.path(lof_context_dir, gene, "cell_line_groups.csv")

  # Screen-specific results remain under the KO output directory.
  gene_dir <- file.path(output_dir, "individual_genes", gene)

  if (!file.exists(group_file)) {
    return(NULL)
  }

  cell_lines <- read_csv(group_file, show_col_types = FALSE)

  ids <- get_screen_ids(
    cell_lines = cell_lines,
    disease = disease,
    mutant_group = mutant_group,
    eligible_ids = eligible_ids
  )

  if (
    length(ids$mutant_id) < min_group_size ||
      length(ids$control_id) < min_group_size
  ) {
    return(NULL)
  }

  intermediate_dir <- file.path(gene_dir, "intermediate_files")
  result_dir <- file.path(gene_dir, "GI_results")
  plot_dir <- file.path(gene_dir, "GI_plots")

  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  if (make_plots) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  context_name <- paste(
    gene,
    mutant_group,
    safe_name(disease),
    sep = "_"
  )

  result_file <- file.path(
    result_dir,
    paste0(
      context_name,
      "_GI_perms_annotated_prioritised_",
      ms_label,
      ".csv"
    )
  )

  # Skip screens that have already been completed.
  if (file.exists(result_file)) {
    message("Skipping existing result: ", basename(result_file))
    return(result_file)
  }

  message(
    "Running ", ms_label, " KO screen: ",
    gene, " | ", mutant_group, " | ", disease,
    " (", length(ids$mutant_id), " mutant; ",
    length(ids$control_id), " control)"
  )

  screen_results <- GRETTA::GI_screen(
    control_id = ids$control_id,
    mutant_id = ids$mutant_id,
    core_num = core_num,
    output_dir = output_dir,
    filename = file.path(
      "individual_genes",
      gene,
      "intermediate_files",
      paste0(context_name, "_GI_", ms_label)
    ),
    data_dir = gretta_data_dir,
    test = FALSE
  )

  perms <- GRETTA::GI_screen_perms(
    control_id = ids$control_id,
    mutant_id = ids$mutant_id,
    n_perm = n_perm,
    core_num = core_num,
    output_dir = output_dir,
    filename = file.path(
      "individual_genes",
      gene,
      "intermediate_files",
      paste0(context_name, "_GI_perms_", ms_label)
    ),
    data_dir = gretta_data_dir
  )

  screen_results_permed <- GRETTA::Append_perm_pval(
    screen_results,
    perms
  )

  intermediate_file <- file.path(
    intermediate_dir,
    paste0(
      context_name,
      "_GI_perms_annotated_",
      ms_label,
      ".csv"
    )
  )

  write_csv(screen_results_permed, intermediate_file)

  screen_results_annotated <- annotate_ko_results(
    screen_results_permed
  )

  write_csv(screen_results_annotated, result_file)

  if (make_plots) {

    plot_data <- screen_results_annotated %>%
      mutate(
        Interaction_score =
          -log10(Pval_perm_adj) * sign(log2FC_by_median)
      )

    genes_to_label <- plot_data %>%
      filter(GI_tier %in% c("Tier_1", "Tier_2", "Tier_3")) %>%
      pull(GeneNames)

    screen_plot <- if (length(genes_to_label) == 0) {
      GRETTA::plot_screen(
        result_df = plot_data,
        label_genes = FALSE
      )
    } else {
      GRETTA::plot_screen(
        result_df = plot_data,
        label_genes = TRUE,
        gene_list = genes_to_label
      )
    }

    ggsave(
      filename = file.path(
        plot_dir,
        paste0(context_name, "_GI_", ms_label, ".pdf")
      ),
      plot = screen_plot,
      width = 8,
      height = 6
    )
  }

  result_file
}


collect_tier123_results <- function(
    screenable_contexts,
    eligible_ids,
    ms_label,
    output_dir) {

  tier123_results <- screenable_contexts %>%
    filter(is_screenable) %>%
    pmap_dfr(function(
        GeneName,
        mut_group,
        disease,
        n_control,
        n_homdel,
        n_thet_del,
        n_mutant,
        is_screenable) {

      result_file <- file.path(
        output_dir,
        "individual_genes",
        GeneName,
        "GI_results",
        paste0(
          GeneName,
          "_",
          mut_group,
          "_",
          safe_name(disease),
          "_GI_perms_annotated_prioritised_",
          ms_label,
          ".csv"
        )
      )

      if (!file.exists(result_file)) {
        return(tibble())
      }

      read_csv(result_file, show_col_types = FALSE) %>%
        filter(GI_tier %in% c("Tier_1", "Tier_2", "Tier_3")) %>%
        mutate(
          Query_gene = GeneName,
          Query_gene_group = mut_group,
          Query_disease_context = disease,
          MS_status = ms_label,
          Total_mutant_line_N = n_mutant,
          Total_control_line_N = n_control,
          Available_HomDel_line_N = n_homdel,
          Available_T_HetDel_line_N = n_thet_del,
          MSI_mutant_line_N = if_else(
            ms_label == "MSI",
            n_mutant,
            0L
          ),
          MSI_control_line_N = if_else(
            ms_label == "MSI",
            n_control,
            0L
          )
        )
    })

  tier123_results
}

# ------------------------------------------------------------------------------
# Run the workflow for MSS and MSI cell lines
# ------------------------------------------------------------------------------

all_tier123_results <- list()

for (ms_label in names(ms_statuses)) {

  inferred_status <- ms_statuses[[ms_label]]

  message("\n", strrep("=", 78))
  message("Processing ", ms_label, " cell lines")
  message(strrep("=", 78))

  eligible_ids <- ms_status_df %>%
    filter(
      Inferred_MS_status == inferred_status,
      DepMap_ID %in% dep$DepMap_ID
    ) %>%
    pull(DepMap_ID) %>%
    unique()

  if (length(eligible_ids) == 0) {
    warning(
      "No cell lines found for microsatellite status: ",
      inferred_status
    )
    next
  }

  # --------------------------------------------------------------------------
  # 1. Identify screenable gene / LOF-group / disease combinations
  # --------------------------------------------------------------------------

  screenable_contexts <- identify_screenable_contexts(
    lof_genes = pog_blof_genes,
    eligible_ids = eligible_ids,
    disease_list = disease_list,
    lof_context_dir = lof_context_dir,
    min_group_size = min_group_size
  )

  screenable_file <- file.path(
    output_dir,
    paste0(
      "KO_screenable_POG_BLOF_contexts_",
      ms_label,
      ".csv"
    )
  )

  write_csv(screenable_contexts, screenable_file)

  message(
    sum(screenable_contexts$is_screenable),
    " screenable ", ms_label, " contexts identified."
  )

  # --------------------------------------------------------------------------
  # 2. Run GRETTA KO screens
  # --------------------------------------------------------------------------

  contexts_to_screen <- screenable_contexts %>%
    filter(is_screenable)

  if (nrow(contexts_to_screen) > 0) {

    walk(
      seq_len(nrow(contexts_to_screen)),
      function(i) {

        context <- contexts_to_screen[i, ]

        run_one_ko_screen(
          gene = context$GeneName,
          mutant_group = context$mut_group,
          disease = context$disease,
          eligible_ids = eligible_ids,
          ms_label = ms_label,
          lof_context_dir = lof_context_dir,
          output_dir = output_dir,
          gretta_data_dir = gretta_data_dir,
          min_group_size = min_group_size,
          n_perm = n_perm,
          core_num = core_num,
          make_plots = make_plots
        )
      }
    )
  }

  # --------------------------------------------------------------------------
  # 3. Collect Tier 1/2/3 results
  # --------------------------------------------------------------------------

  tier123_results <- collect_tier123_results(
    screenable_contexts = screenable_contexts,
    eligible_ids = eligible_ids,
    ms_label = ms_label,
    output_dir = output_dir
  )

  tier123_file <- file.path(
    output_dir,
    paste0(
      "KO_POG_BLOF_",
      ms_label,
      "_screen_results_Tier123.csv"
    )
  )

  write_csv(tier123_results, tier123_file)

  all_tier123_results[[ms_label]] <- tier123_results
}

# ------------------------------------------------------------------------------
# Combine MSS and MSI Tier 1/2/3 results
# ------------------------------------------------------------------------------

combined_tier123_results <- bind_rows(all_tier123_results)

write_csv(
  combined_tier123_results,
  file.path(
    output_dir,
    "KO_POG_BLOF_MSS_MSI_screen_results_Tier123.csv"
  )
)

message("\nKO screening workflow complete.")
