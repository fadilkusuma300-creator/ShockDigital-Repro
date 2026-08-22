args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
dir0 <- if (length(self0)) dirname(normalizePath(self0[1], mustWork = TRUE)) else "scripts"
source(file.path(dir0, "_bootstrap.R"))
source(file.path(root, "R", "08_balance_tables.R"))
source(file.path(root, "R", "11_figures_tables.R"))

cfg <- load_analysis_config(root)
ensure_project_dirs(cfg, root)
windows <- readRDS(file.path(root, cfg$paths$derived_dir, "strict_three_wave_sample.rds"))
bal <- data.table::fread(file.path(root, cfg$paths$results_dir, "diagnostics", "balance_all_covariates.csv"))

# Figure 1 uses the Python framework generator when available, with an R fallback.
fig1 <- file.path(root, cfg$paths$results_dir, "figures", "Fig1.pdf")
py <- Sys.which("python"); if (!nzchar(py)) py <- Sys.which("python3")
fig1_script <- file.path(root, "python", "fig1_framework.py")
fig_cfg <- file.path(root, "config", "figures.yml")
status <- 1L
if (nzchar(py) && file.exists(fig1_script) && file.exists(fig_cfg)) {
  status <- system2(py, c(fig1_script, "--config", fig_cfg, "--output", fig1))
}
if (!identical(status, 0L) || !file.exists(fig1)) save_fig1_design(fig1)

# Figures 2-4 are rendered from the saved numerical outputs.
curve_path <- file.path(root, cfg$paths$results_dir, "tables", "primary_sales_curve.csv")
thr_path <- file.path(root, cfg$paths$results_dir, "tables", "sustained_support_threshold.csv")
sev_path <- file.path(root, cfg$paths$results_dir, "tables", "severe_shock_outcomes.csv")
sec_path <- file.path(root, cfg$paths$results_dir, "tables", "secondary_effect_curves.csv")
if (file.exists(curve_path) && file.exists(thr_path)) {
  curve <- data.table::fread(curve_path)
  threshold <- data.table::fread(thr_path)$sustained_support_threshold[1]
  save_fig2_primary_curve(curve, threshold, file.path(root, cfg$paths$results_dir, "figures", "Fig2.pdf"))
}
if (file.exists(sev_path)) {
  severe <- data.table::fread(sev_path)
  save_fig3_severe_forest(severe, file.path(root, cfg$paths$results_dir, "figures", "Fig3.pdf"))
}
if (file.exists(sec_path)) {
  curves <- data.table::fread(sec_path)
  save_fig4_secondary_curves(curves, file.path(root, cfg$paths$results_dir, "figures", "Fig4.pdf"))
}

panel_b <- sample_panel_b(windows, bal)
write_csv_safely(panel_b, file.path(root, cfg$paths$results_dir, "tables", "sample_balance_panel_b.csv"))
