args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
dir0 <- if (length(self0)) dirname(normalizePath(self0[1], mustWork = TRUE)) else "scripts"
source(file.path(dir0, "_bootstrap.R"))
source(file.path(root, "R", "10_validation.R"))

cfg <- load_analysis_config(root)
primary <- readRDS(file.path(root, cfg$paths$derived_dir, "primary_sales_dr.rds"))
reps_env <- suppressWarnings(as.integer(Sys.getenv("VALIDATION_REPS", "")))
reps <- if (is.finite(reps_env) && reps_env > 0) reps_env else cfg$validation$repetitions
methods_env <- Sys.getenv("VALIDATION_METHODS", "")
methods <- if (nzchar(methods_env)) strsplit(methods_env, ",", fixed = TRUE)[[1]] else
  c("grf", "causalpfn", "model_averaged_dml", "conditional_dr")

raw <- list()
for (mechanism in c("smooth_transition", "sharp_jump", "linear_change")) {
  for (r in seq_len(reps)) {
    log_message("Validation ", mechanism, " repetition ", r, "/", reps, root = root)
    raw[[length(raw) + 1L]] <- run_validation_repetition(primary, mechanism, r, cfg, methods)
  }
}
raw <- data.table::rbindlist(raw)
summary <- summarize_validation(raw)
write_csv_safely(raw, file.path(root, cfg$paths$results_dir, "validation", "validation_raw.csv"))
write_csv_safely(summary, file.path(root, cfg$paths$results_dir, "tables", "validation_summary.csv"))

component <- data.table::rbindlist(lapply(seq_len(reps), function(r) {
  log_message("Component analysis repetition ", r, "/", reps, root = root)
  run_component_repetition(primary, r, cfg)
}))
component_summary <- summarize_validation(component)
write_csv_safely(component, file.path(root, cfg$paths$results_dir, "validation", "component_raw.csv"))
write_csv_safely(component_summary, file.path(root, cfg$paths$results_dir, "tables", "component_summary.csv"))
