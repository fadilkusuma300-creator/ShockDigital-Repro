args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
dir0 <- if (length(self0)) dirname(normalizePath(self0[1], mustWork = TRUE)) else "scripts"
source(file.path(dir0, "_bootstrap.R"))
source(file.path(root, "R", "09_robustness.R"))

cfg <- load_analysis_config(root)
windows <- readRDS(file.path(root, cfg$paths$derived_dir, "strict_three_wave_sample.rds"))
primary <- readRDS(file.path(root, cfg$paths$derived_dir, "primary_sales_dr.rds"))

log_message("Running robustness analyses", root = root)
rob <- run_sales_robustness(windows, primary, cfg)
write_csv_safely(rob$summary, file.path(root, cfg$paths$results_dir, "tables", "robustness_summary.csv"))
write_csv_safely(rob$drop_one_economy, file.path(root, cfg$paths$results_dir, "diagnostics", "leave_one_economy_out.csv"))
save_rds_safely(rob, file.path(root, cfg$paths$derived_dir, "robustness_results.rds"))
