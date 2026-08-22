# Install the R packages required by the analysis pipeline.

pkgs <- c(
  "data.table", "haven", "yaml", "ranger", "glmnet",
  "xgboost", "quadprog", "sandwich", "ggplot2", "grf"
)

repos <- getOption("repos")
if (is.null(repos) || identical(repos[["CRAN"]], "@CRAN@")) {
  repos[["CRAN"]] <- "https://cloud.r-project.org"
}

missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) {
  install.packages(missing, repos = repos, dependencies = TRUE)
}

message("R dependencies are installed.")
