args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
dir0 <- if (length(self0)) dirname(normalizePath(self0[1], mustWork = TRUE)) else "scripts"
source(file.path(dir0, "_bootstrap.R"))
source(file.path(root, "R", "03_preprocess.R"))
source(file.path(root, "R", "04_nuisance_models.R"))
source(file.path(root, "R", "05_dr_estimator.R"))
source(file.path(root, "R", "06_continuous_curve.R"))
source(file.path(root, "R", "07_inference.R"))
source(file.path(root, "R", "11_figures_tables.R"))

cfg <- load_analysis_config(root)
windows <- readRDS(file.path(root, cfg$paths$derived_dir, "strict_three_wave_sample.rds"))
primary <- readRDS(file.path(root, cfg$paths$derived_dir, "primary_sales_dr.rds"))

outcomes <- c("y_sales", "y_continued_operation", "y_employment_retention", "y_input_recovery")
fits <- list(y_sales = primary)
for (j in seq_along(outcomes[-1])) {
  outcome <- outcomes[-1][j]
  log_message("Estimating secondary outcome: ", outcome, root = root)
  fits[[outcome]] <- crossfit_dr_pseudo_outcome(windows, outcome, cfg, seed_offset = 50L + j)
}

severe <- severe_outcome_table(fits, cfg)
write_csv_safely(severe, file.path(root, cfg$paths$results_dir, "tables", "severe_shock_outcomes.csv"))
save_fig3_severe_forest(severe, file.path(root, cfg$paths$results_dir, "figures", "Fig3.pdf"))

shock_grid <- seq(cfg$sample$primary_curve_min, cfg$sample$primary_curve_max,
                  by = cfg$sample$shock_grid_step)
secondary_curves <- list()
for (j in seq_along(outcomes[-1])) {
  outcome <- outcomes[-1][j]
  m <- fit_continuous_effect_curve(fits[[outcome]]$data, cfg, seed_offset = 100L + j)
  ctab <- pointwise_curve(m, fits[[outcome]]$data, shock_grid, cfg)
  # Binary recovery outcomes are reported in percentage points.
  ctab[, `:=`(estimate = 100 * estimate, se = 100 * se,
              lower = 100 * lower, upper = 100 * upper)]
  ctab[, outcome := outcome]
  secondary_curves[[outcome]] <- ctab
  save_rds_safely(m, file.path(root, cfg$paths$derived_dir, paste0(outcome, "_curve_model.rds")))
}
curves <- data.table::rbindlist(secondary_curves)
write_csv_safely(curves, file.path(root, cfg$paths$results_dir, "tables", "secondary_effect_curves.csv"))
save_fig4_secondary_curves(curves, file.path(root, cfg$paths$results_dir, "figures", "Fig4.pdf"))
save_rds_safely(fits, file.path(root, cfg$paths$derived_dir, "all_outcome_dr_fits.rds"))
