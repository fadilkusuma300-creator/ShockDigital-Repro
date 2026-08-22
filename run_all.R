# End-to-end analysis pipeline.
args0 <- commandArgs(trailingOnly = FALSE)
self0 <- sub("^--file=", "", args0[grep("^--file=", args0)])
if (length(self0)) {
  root0 <- normalizePath(dirname(self0[1]), winslash = "/", mustWork = TRUE)
  setwd(root0)
}

rscript <- Sys.which("Rscript")
if (!nzchar(rscript)) stop("Rscript executable not found")

stages <- sprintf("scripts/%02d_%s.R", 1:6, c(
  "prepare_data", "primary_sales", "secondary_outcomes",
  "robustness", "validation", "figures_tables"
))

for (stage in stages) {
  message("Running ", stage)
  status <- system2(rscript, stage)
  if (!identical(status, 0L)) stop("Analysis stage failed: ", stage)
}
