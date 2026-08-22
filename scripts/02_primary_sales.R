args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
dir0 <- if (length(self0)) dirname(normalizePath(self0[1], mustWork = TRUE)) else "scripts"
source(file.path(dir0, "_bootstrap.R"))
source(file.path(root, "R", "03_preprocess.R"))
source(file.path(root, "R", "04_nuisance_models.R"))
source(file.path(root, "R", "05_dr_estimator.R"))
source(file.path(root, "R", "06_continuous_curve.R"))
source(file.path(root, "R", "07_inference.R"))
source(file.path(root, "R", "08_balance_tables.R"))
source(file.path(root, "R", "11_figures_tables.R"))

cfg <- load_analysis_config(root)
ensure_project_dirs(cfg, root)
windows <- readRDS(file.path(root, cfg$paths$derived_dir, "strict_three_wave_sample.rds"))

log_message("Estimating primary sales model", root = root)
primary <- crossfit_dr_pseudo_outcome(
  windows, "y_sales", cfg,
  support = c(cfg$sample$propensity_lower, cfg$sample$propensity_upper)
)
ate <- estimate_ate_table(primary, cfg)
bal <- balance_table(primary$data)
bal_summary <- balance_summary(bal)

shock_grid <- seq(
  cfg$sample$primary_curve_min,
  cfg$sample$primary_curve_max,
  by = cfg$sample$shock_grid_step
)
curve_model <- fit_continuous_effect_curve(primary$data, cfg)
curve <- wild_cluster_simultaneous_band(curve_model, primary$data, shock_grid, cfg)
threshold <- sustained_support_point(curve, cfg)

save_rds_safely(primary, file.path(root, cfg$paths$derived_dir, "primary_sales_dr.rds"))
save_rds_safely(curve_model, file.path(root, cfg$paths$derived_dir, "primary_sales_curve_model.rds"))
write_csv_safely(ate, file.path(root, cfg$paths$results_dir, "tables", "primary_ate.csv"))
write_csv_safely(bal, file.path(root, cfg$paths$results_dir, "diagnostics", "balance_all_covariates.csv"))
write_csv_safely(bal_summary, file.path(root, cfg$paths$results_dir, "diagnostics", "balance_summary.csv"))
write_csv_safely(primary$ensemble_weights, file.path(root, cfg$paths$results_dir, "diagnostics", "outcome_ensemble_weights_sales.csv"))
write_csv_safely(curve_model$tuning_grid, file.path(root, cfg$paths$results_dir, "diagnostics", "second_stage_tuning_sales.csv"))
write_csv_safely(curve, file.path(root, cfg$paths$results_dir, "tables", "primary_sales_curve.csv"))
write_csv_safely(data.table::data.table(sustained_support_threshold = threshold),
                 file.path(root, cfg$paths$results_dir, "tables", "sustained_support_threshold.csv"))

save_fig2_primary_curve(
  curve, threshold,
  file.path(root, cfg$paths$results_dir, "figures", "Fig2.pdf")
)
log_message("Primary sales model complete; sustained-support point = ", threshold, root = root)
