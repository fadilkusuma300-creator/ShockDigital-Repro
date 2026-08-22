args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
dir0 <- if (length(self0)) dirname(normalizePath(self0[1], mustWork = TRUE)) else "scripts"
source(file.path(dir0, "_bootstrap.R"))
source(file.path(root, "R", "01_io_harmonize.R"))
source(file.path(root, "R", "02_sample_windows.R"))
source(file.path(root, "R", "03_preprocess.R"))

cfg <- load_analysis_config(root)
vmap <- load_variable_map(root)
ensure_project_dirs(cfg, root)
set_project_seed(cfg)
log_message("Starting data harmonization", root = root)

baseline <- harmonize_baseline(root, cfg, vmap)
followup <- harmonize_followup(root, cfg, vmap)
followup <- prepare_survey_weights(followup)
followup <- add_relative_interview_month(followup)

# Keep an all-size matched panel for the sample-flow count, then apply the
# 5-249 baseline-employment definition for all analytical work.
matched_all <- merge_baseline_followup(baseline, followup, cfg, apply_sme_filter = FALSE)
merged_sme <- merge_baseline_followup(baseline, followup, cfg, apply_sme_filter = TRUE)
windows <- construct_adjacent_windows(merged_sme, cfg)
windows <- prepare_raw_analysis_columns(windows)

flow <- sample_flow_table(baseline, followup, matched_all, merged_sme, windows, cfg)
write_csv_safely(flow, file.path(root, cfg$paths$results_dir, "tables", "sample_flow.csv"))
detail_flow <- detailed_sample_flow_table(merged_sme, windows)
write_csv_safely(detail_flow, file.path(root, cfg$paths$results_dir, "diagnostics", "sample_flow_detailed.csv"))
save_rds_safely(windows, file.path(root, cfg$paths$derived_dir, "strict_three_wave_sample.rds"))
save_rds_safely(followup, file.path(root, cfg$paths$derived_dir, "harmonized_followup.rds"))
save_rds_safely(baseline, file.path(root, cfg$paths$derived_dir, "harmonized_baseline.rds"))

# Baseline finance and linkage summaries.
credit_summary <- baseline[, .(
  n = .N,
  classified = sum(!is.na(credit_constraints)),
  constrained = sum(credit_constraints == 1L, na.rm = TRUE),
  unconstrained = sum(credit_constraints == 0L, na.rm = TRUE),
  missing = sum(is.na(credit_constraints))
), by = .(source_file, economy, credit_constraint_source)][order(economy, source_file)]
write_csv_safely(credit_summary, file.path(root, cfg$paths$results_dir, "diagnostics", "credit_constraint_summary.csv"))

baseline_fields <- c(
  "baseline_employment", "firm_age", "manager_experience", "website",
  "credit_constraints", "innovation", "foreign_owned", "exporter", "industry"
)
baseline_fields <- intersect(baseline_fields, names(matched_all))
complete_flag <- rep(TRUE, nrow(matched_all))
for (nm in baseline_fields) {
  z <- matched_all[[nm]]
  if (is.numeric(z) || is.integer(z)) complete_flag <- complete_flag & is.finite(safe_numeric(z))
  else complete_flag <- complete_flag & !is.na(z) & nzchar(trimws(as.character(z)))
}
baseline_summary <- data.table::data.table(
  metric = c("matched_to_baseline_unique_firms", "complete_core_baseline_before_imputation_unique_firms"),
  observed = c(
    data.table::uniqueN(matched_all, by = c("economy", "firm_id")),
    data.table::uniqueN(matched_all[complete_flag], by = c("economy", "firm_id"))
  ),
  note = c(
    "Linked to baseline; individual covariates may still contain missing values.",
    paste0("Complete on: ", paste(baseline_fields, collapse = ", "), ".")
  )
)
write_csv_safely(baseline_summary, file.path(root, cfg$paths$results_dir, "diagnostics", "baseline_completeness.csv"))

log_message("Prepared strict sample with ", nrow(windows), " firms across ",
            data.table::uniqueN(windows$economy), " economies", root = root)
