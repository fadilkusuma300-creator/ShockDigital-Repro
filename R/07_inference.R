# Cluster-robust pointwise inference and wild-cluster simultaneous confidence band.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "06_continuous_curve.R"))

curve_influence_weights <- function(model, target_data, shock_grid) {
  cmat <- standardization_rows(model, target_data, shock_grid)
  # theta = A^-1 X' W y; tau(s)=c(s)'theta.
  # Thus the observation-level linear weight is W X A^-1 c(s).
  h <- model$design %*% model$normal_inv %*% t(cmat)
  h <- h * model$weights
  list(cmat = cmat, h = h)
}

cluster_curve_se <- function(model, target_data, shock_grid) {
  ih <- curve_influence_weights(model, target_data, shock_grid)
  h <- ih$h
  r <- model$residuals
  clusters <- as.character(model$data$economy)
  ug <- unique(clusters)
  g <- length(ug)
  if (g < 2L) stop("At least two economy clusters are required for inference")

  scores <- matrix(0, g, length(shock_grid))
  for (j in seq_along(ug)) {
    idx <- which(clusters == ug[j])
    scores[j, ] <- colSums(h[idx, , drop = FALSE] * r[idx])
  }
  correction <- g / (g - 1)
  se <- sqrt(correction * colSums(scores^2))
  list(se = se, cluster_scores = scores, influence = ih)
}

wild_cluster_simultaneous_band <- function(model, target_data, shock_grid, cfg,
                                           seed_offset = 0L) {
  base <- as.numeric(standardization_rows(model, target_data, shock_grid) %*% model$beta)
  inf <- cluster_curve_se(model, target_data, shock_grid)
  se <- inf$se
  clusters <- as.character(model$data$economy)
  ug <- unique(clusters)
  gindex <- match(clusters, ug)
  b <- cfg$inference$wild_cluster_replications
  set.seed(cfg$seed + 7000L + seed_offset)

  cmat <- inf$influence$cmat
  x <- model$design
  w <- model$weights
  ainv <- model$normal_inv
  fitted <- model$fitted
  resid <- model$residuals
  max_t <- numeric(b)

  for (r in seq_len(b)) {
    mult <- sample(c(-1, 1), length(ug), replace = TRUE)
    y_star <- fitted + resid * mult[gindex]
    rhs <- crossprod(x, w * y_star)
    beta_star <- as.numeric(ainv %*% rhs)
    curve_star <- as.numeric(cmat %*% beta_star)
    tstat <- abs((curve_star - base) / pmax(se, 1e-12))
    max_t[r] <- max(tstat[is.finite(tstat)], na.rm = TRUE)
  }

  alpha <- 1 - cfg$inference$pointwise_level
  crit <- as.numeric(stats::quantile(max_t, 1 - alpha, na.rm = TRUE, type = 8))
  z <- stats::qnorm(1 - alpha / 2)
  data.table::data.table(
    shock = shock_grid,
    estimate = base,
    se = se,
    pointwise_lower = base - z * se,
    pointwise_upper = base + z * se,
    simultaneous_lower = base - crit * se,
    simultaneous_upper = base + crit * se,
    simultaneous_critical_value = crit
  )
}

sustained_support_point <- function(curve, cfg) {
  k <- as.integer(cfg$inference$sustained_support_adjacent_points)
  lower <- curve$simultaneous_lower
  for (i in seq_len(nrow(curve))) {
    j <- i + k
    if (j <= nrow(curve) && all(is.finite(lower[i:j])) && all(lower[i:j] > 0)) {
      return(curve$shock[i])
    }
  }
  NA_real_
}

pointwise_curve <- function(model, target_data, shock_grid, cfg) {
  est <- as.numeric(standardization_rows(model, target_data, shock_grid) %*% model$beta)
  se <- cluster_curve_se(model, target_data, shock_grid)$se
  z <- stats::qnorm(1 - (1 - cfg$inference$pointwise_level) / 2)
  data.table::data.table(
    shock = shock_grid, estimate = est, se = se,
    lower = est - z * se, upper = est + z * se
  )
}
