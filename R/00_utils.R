# Shared utility functions.
# All functions are deterministic when the project seed is fixed.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

project_root <- function() {
  # Resolve the project directory whether a script is run from the root,
  # scripts/, R/, or an IDE working directory.
  candidates <- c(".", "..", "../..")
  for (p in candidates) {
    if (file.exists(file.path(p, "config", "analysis.yml"))) {
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Cannot locate project root containing config/analysis.yml")
}

load_analysis_config <- function(root = project_root()) {
  yaml::read_yaml(file.path(root, "config", "analysis.yml"))
}

load_variable_map <- function(root = project_root()) {
  yaml::read_yaml(file.path(root, "config", "variable_map.yml"))
}

load_figure_config <- function(root = project_root()) {
  yaml::read_yaml(file.path(root, "config", "figures.yml"))
}


ensure_project_dirs <- function(cfg, root = project_root()) {
  dirs <- c(
    file.path(root, cfg$paths$derived_dir),
    file.path(root, cfg$paths$results_dir),
    file.path(root, cfg$paths$results_dir, "figures"),
    file.path(root, cfg$paths$results_dir, "tables"),
    file.path(root, cfg$paths$results_dir, "diagnostics"),
    file.path(root, cfg$paths$results_dir, "validation"),
    file.path(root, "logs")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

set_project_seed <- function(cfg, offset = 0L) {
  seed <- as.integer(cfg$seed) + as.integer(offset)
  set.seed(seed)
  invisible(seed)
}

clip <- function(x, lower, upper) {
  pmin(pmax(x, lower), upper)
}

normalize_weights <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w <= 0] <- NA_real_
  if (all(is.na(w))) return(rep(1, length(w)))
  w[is.na(w)] <- stats::median(w, na.rm = TRUE)
  w / mean(w)
}

weighted_mean <- function(x, w, na.rm = TRUE) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2L) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  mu <- sum(w * x) / sum(w)
  denom <- sum(w) - sum(w^2) / sum(w)
  if (denom <= 0) return(NA_real_)
  sum(w * (x - mu)^2) / denom
}

weighted_sd <- function(x, w) sqrt(weighted_var(x, w))

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(rep(NA_real_, length(probs)))
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) {
    idx <- which(cw >= p)[1]
    x[idx]
  }, numeric(1))
}

weighted_mode <- function(x, w) {
  ok <- !is.na(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA)
  tab <- tapply(w[ok], x[ok], sum)
  names(tab)[which.max(tab)]
}

winsorize_weight <- function(w, probs = c(0.01, 0.99)) {
  q <- weighted_quantile(w, rep(1, length(w)), probs)
  clip(w, q[1], q[2])
}

safe_numeric <- function(x) {
  if (inherits(x, "haven_labelled")) x <- haven::zap_labels(x)
  suppressWarnings(as.numeric(x))
}

safe_character <- function(x) {
  if (inherits(x, "haven_labelled")) x <- haven::as_factor(x, levels = "values")
  as.character(x)
}

missing_to_na <- function(x, missing_codes) {
  if (is.numeric(x) || is.integer(x) || inherits(x, "haven_labelled")) {
    z <- safe_numeric(x)
    z[z %in% unlist(missing_codes)] <- NA_real_
    return(z)
  }
  x
}

make_stratified_folds <- function(economy, treatment, k = 5L, seed = 1L) {
  # Stratification is performed on economy x treatment cells whenever enough
  # observations are available. Small cells are distributed cyclically so no
  # observation is discarded merely because an economy-treatment cell is tiny.
  set.seed(seed)
  n <- length(treatment)
  folds <- integer(n)
  strata <- interaction(as.character(economy), as.character(treatment), drop = TRUE)
  for (lev in levels(strata)) {
    idx <- which(strata == lev)
    idx <- sample(idx, length(idx), replace = FALSE)
    folds[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  # If any fold is empty because the data are extremely small, use a global
  # cyclic fallback while preserving randomness.
  if (any(tabulate(folds, nbins = k) == 0L)) {
    idx <- sample(seq_len(n), n, replace = FALSE)
    folds[idx] <- rep(seq_len(k), length.out = n)
  }
  folds
}

make_group_folds <- function(group, k = 3L, seed = 1L) {
  # Keeps entire economies together in the same fold. This is useful for
  # second-stage hyperparameter selection so cluster structure is respected.
  set.seed(seed)
  ug <- sample(unique(as.character(group)))
  gfold <- setNames(rep(seq_len(k), length.out = length(ug)), ug)
  unname(gfold[as.character(group)])
}

weighted_smd_numeric <- function(x, a, w) {
  i1 <- which(a == 1 & is.finite(x) & is.finite(w))
  i0 <- which(a == 0 & is.finite(x) & is.finite(w))
  if (length(i1) < 2L || length(i0) < 2L) return(NA_real_)
  m1 <- weighted_mean(x[i1], w[i1])
  m0 <- weighted_mean(x[i0], w[i0])
  v1 <- weighted_var(x[i1], w[i1])
  v0 <- weighted_var(x[i0], w[i0])
  den <- sqrt((v1 + v0) / 2)
  if (!is.finite(den) || den == 0) return(0)
  (m1 - m0) / den
}

weighted_smd_binary <- function(x, a, w) {
  x <- as.numeric(x)
  weighted_smd_numeric(x, a, w)
}

cluster_mean_se <- function(signal, weights, cluster) {
  # Influence-function standard error for a weighted mean clustered by economy.
  ok <- is.finite(signal) & is.finite(weights) & weights > 0 & !is.na(cluster)
  signal <- signal[ok]
  weights <- weights[ok]
  cluster <- as.character(cluster[ok])
  mu <- sum(weights * signal) / sum(weights)
  infl <- weights * (signal - mu) / sum(weights)
  cg <- tapply(infl, cluster, sum)
  g <- length(cg)
  n <- length(signal)
  if (g < 2L) return(list(estimate = mu, se = NA_real_, clusters = g, n = n))
  correction <- (g / (g - 1)) * ((n - 1) / max(1, n - 1))
  se <- sqrt(correction * sum(cg^2))
  list(estimate = mu, se = se, clusters = g, n = n)
}

ci_from_estimate <- function(est, se, level = 0.95) {
  alpha <- 1 - level
  z <- stats::qnorm(1 - alpha / 2)
  c(lower = est - z * se, upper = est + z * se)
}

log_message <- function(..., root = project_root()) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  cat(txt, "\n")
  cat(txt, "\n", file = file.path(root, "logs", "analysis.log"), append = TRUE)
}

write_csv_safely <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(data.table::as.data.table(x), path)
  invisible(path)
}

save_rds_safely <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
  invisible(path)
}
