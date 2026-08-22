# Resolve the project root for command-line and IDE execution.
args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_full[grep("^--file=", args_full)])
if (length(file_arg)) {
  root <- normalizePath(file.path(dirname(file_arg[1]), ".."), winslash = "/", mustWork = TRUE)
} else {
  candidates <- c(".", "..")
  root <- NULL
  for (p in candidates) {
    if (file.exists(file.path(p, "config", "analysis.yml"))) {
      root <- normalizePath(p, winslash = "/", mustWork = TRUE)
      break
    }
  }
  if (is.null(root)) stop("Cannot locate project root")
}
setwd(root)
source(file.path(root, "R", "00_utils.R"))
