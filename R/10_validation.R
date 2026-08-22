# Semi-synthetic validation and component analysis.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "03_preprocess.R"))
source(file.path(project_root(), "R", "05_dr_estimator.R"))
source(file.path(project_root(), "R", "06_continuous_curve.R"))
source(file.path(project_root(), "R", "07_inference.R"))

calibrate_propensity_share <- function(p, w = NULL, target = 0.24) {
  p <- clip(p, 1e-5, 1 - 1e-5)
  f <- function(delta) {
    q <- stats::plogis(stats::qlogis(p) + delta)
    mean(q) - target
  }
  delta <- stats::uniroot(f, interval = c(-12, 12))$root
  stats::plogis(stats::qlogis(p) + delta)
}

shock_quintile <- function(s, w, groups = 5L) {
  q <- weighted_quantile(s, w, seq(0, 1, length.out = groups + 1L))
  q[1] <- -Inf
  q[length(q)] <- Inf
  # Duplicate cutpoints are separated by a tiny amount only for bin creation.
  for (j in 2:(length(q) - 1L)) {
    if (q[j] <= q[j - 1L]) q[j] <- q[j - 1L] + 1e-8
  }
  as.integer(cut(s, breaks = q, include.lowest = TRUE, labels = FALSE))
}

sample_control_residuals <- function(base, seed, groups = 5L) {
  set.seed(seed)
  x <- data.table::copy(base)
  x[, shock_q := shock_quintile(shock_sales_decline, analysis_weight, groups)]
  x[, control_residual := y_sales - m0_hat]
  pools <- x[treatment == 0L & is.finite(control_residual),
             .(pool = list(control_residual)), by = .(economy, shock_q)]
  key <- paste(x$economy, x$shock_q, sep = "||")
  lookup <- setNames(pools$pool, paste(pools$economy, pools$shock_q, sep = "||"))
  missing_cells <- setdiff(unique(key), names(lookup))
  if (length(missing_cells)) {
    stop("Validation residual pool missing for economy x shock-quintile cells: ",
         paste(missing_cells, collapse = ", "))
  }
  vapply(seq_len(nrow(x)), function(i) {
    sample(lookup[[key[i]]], 1L, replace = TRUE)
  }, numeric(1))
}

true_effect_function <- function(s, mechanism, parameters, shock_range = c(0, 80)) {
  if (mechanism == "smooth_transition") {
    # The paper gives center, transition width, and plateau magnitude but not a
    # closed-form transition function. Here width is defined as the 10%-to-90%
    # logistic transition span, which keeps the maximum slope at the stated center.
    center <- parameters$center
    width <- parameters$width
    plateau <- parameters$plateau
    scale <- width / (2 * log(9))
    return(plateau * stats::plogis((s - center) / scale))
  }
  if (mechanism == "sharp_jump") {
    return(parameters$effect * as.numeric(s >= parameters$location))
  }
  if (mechanism == "linear_change") {
    z <- clip((s - shock_range[1]) / diff(shock_range), 0, 1)
    return(parameters$total_effect * z)
  }
  stop("Unknown mechanism: ", mechanism)
}

draw_dgp_parameters <- function(mechanism, cfg) {
  u <- function(r) stats::runif(1, r[1], r[2])
  if (mechanism == "smooth_transition") {
    z <- cfg$validation$smooth_transition
    return(list(center = u(z$center_range), width = u(z$width_range), plateau = u(z$plateau_range)))
  }
  if (mechanism == "sharp_jump") {
    z <- cfg$validation$sharp_jump
    return(list(location = u(z$location_range), effect = u(z$effect_range)))
  }
  if (mechanism == "linear_change") {
    z <- cfg$validation$linear_change
    return(list(total_effect = u(z$total_effect_range)))
  }
  stop("Unknown mechanism: ", mechanism)
}

generate_semisynthetic_data <- function(primary_fit, mechanism, cfg, seed) {
  set.seed(seed)
  x <- data.table::copy(primary_fit$data)
  p <- calibrate_propensity_share(
    x$propensity_raw, x$analysis_weight, cfg$validation$treatment_share_target
  )
  a_sim <- stats::rbinom(nrow(x), size = 1, prob = p)
  eps <- sample_control_residuals(
    x, seed + 1L, cfg$validation$residual_shock_quintiles
  )
  pars <- draw_dgp_parameters(mechanism, cfg)
  tau <- true_effect_function(
    x$shock_sales_decline, mechanism, pars,
    c(cfg$sample$primary_curve_min, cfg$sample$primary_curve_max)
  )
  y0 <- x$m0_hat + eps
  y <- y0 + a_sim * tau
  x[, `:=`(
    treatment = a_sim,
    y_sim = y,
    y0_true = y0,
    tau_true_observed_s = tau,
    simulated_propensity = p
  )]
  attr(x, "dgp_parameters") <- pars
  x
}

curve_metrics <- function(curve, true_curve, mechanism, parameters) {
  m <- merge(curve, true_curve, by = "shock", all = FALSE)
  rng <- diff(range(m$true_effect, na.rm = TRUE))
  rmse <- sqrt(mean((m$estimate - m$true_effect)^2, na.rm = TRUE))
  standardized <- if (rng > 0) rmse / rng else NA_real_
  lower_col <- intersect(c("lower", "pointwise_lower"), names(m))[1]
  upper_col <- intersect(c("upper", "pointwise_upper"), names(m))[1]
  coverage <- NA_real_
  if (!is.na(lower_col) && !is.na(upper_col)) {
    coverage <- mean(
      m$true_effect >= m[[lower_col]] & m$true_effect <= m[[upper_col]],
      na.rm = TRUE
    )
  }

  change_error <- NA_real_
  if (mechanism %in% c("smooth_transition", "sharp_jump") && nrow(m) > 2L) {
    slope <- diff(m$estimate) / diff(m$shock)
    estimated_location <- m$shock[-nrow(m)][which.max(abs(slope))]
    true_location <- if (mechanism == "smooth_transition") parameters$center else parameters$location
    change_error <- abs(estimated_location - true_location)
  }
  data.table::data.table(
    standardized_rmse = standardized,
    change_location_error = change_error,
    pointwise_coverage = coverage
  )
}

fit_primary_validation_method <- function(sim, cfg, shock_grid, seed_offset) {
  dr <- crossfit_dr_pseudo_outcome(sim, "y_sim", cfg, seed_offset = seed_offset)
  model <- fit_continuous_effect_curve(dr$data, cfg, seed_offset = seed_offset)
  pointwise_curve(model, dr$data, shock_grid, cfg)
}

fit_no_orthogonalization <- function(sim, cfg, shock_grid, seed_offset) {
  dr <- crossfit_dr_pseudo_outcome(sim, "y_sim", cfg, seed_offset = seed_offset)
  x <- data.table::copy(dr$data)
  x[, dr_pseudo_outcome := m1_hat - m0_hat]
  model <- fit_continuous_effect_curve(x, cfg, seed_offset = seed_offset)
  pointwise_curve(model, x, shock_grid, cfg)
}

fit_unstandardized_spline <- function(sim, cfg, shock_grid, seed_offset) {
  dr <- crossfit_dr_pseudo_outcome(sim, "y_sim", cfg, seed_offset = seed_offset)
  x <- dr$data
  s <- x$shock_sales_decline
  w <- normalize_weights(x$analysis_weight)
  boundary <- range(s)
  knots <- weighted_quantile(s, w, unlist(cfg$second_stage$spline_knot_probs))
  b <- make_spline_basis(s, knots, boundary, cfg$second_stage$spline_degree)
  bg <- make_spline_basis(shock_grid, knots, boundary, cfg$second_stage$spline_degree)
  lambdas <- unlist(cfg$second_stage$lambda_spline_grid)
  folds <- make_group_folds(x$economy, cfg$second_stage$inner_folds, cfg$seed + seed_offset)

  score <- vapply(lambdas, function(lambda) {
    fold_loss <- numeric(max(folds))
    for (f in seq_len(max(folds))) {
      tr <- which(folds != f); va <- which(folds == f)
      p <- lambda * second_difference_penalty(ncol(b))
      a <- crossprod(b[tr, , drop = FALSE], w[tr] * b[tr, , drop = FALSE]) + p + diag(1e-10, ncol(b))
      beta <- solve(a, crossprod(b[tr, , drop = FALSE], w[tr] * x$dr_pseudo_outcome[tr]))
      pred <- as.numeric(b[va, , drop = FALSE] %*% beta)
      fold_loss[f] <- weighted_mean((x$dr_pseudo_outcome[va] - pred)^2, w[va])
    }
    mean(fold_loss)
  }, numeric(1))
  lambda <- lambdas[which.min(score)]
  p <- lambda * second_difference_penalty(ncol(b))
  normal <- crossprod(b, w * b) + p + diag(1e-10, ncol(b))
  beta <- solve(normal, crossprod(b, w * x$dr_pseudo_outcome))
  est <- as.numeric(bg %*% beta)

  # Cluster sandwich for the penalized linear smoother.
  residual <- x$dr_pseudo_outcome - as.numeric(b %*% beta)
  ainv <- solve(normal)
  h <- (b %*% ainv %*% t(bg)) * w
  clusters <- as.character(x$economy)
  ug <- unique(clusters)
  scores <- matrix(0, length(ug), length(shock_grid))
  for (j in seq_along(ug)) {
    idx <- which(clusters == ug[j])
    scores[j, ] <- colSums(h[idx, , drop = FALSE] * residual[idx])
  }
  se <- sqrt(length(ug) / (length(ug) - 1) * colSums(scores^2))
  z <- stats::qnorm(1 - (1 - cfg$inference$pointwise_level) / 2)
  data.table::data.table(shock = shock_grid, estimate = est,
                         lower = est - z * se, upper = est + z * se)
}

fit_model_averaged_dml_curve <- function(sim, cfg, shock_grid, seed_offset) {
  # Residualized DML curve using the same cross-fitted nuisance learners, followed
  # by a spline coefficient model for the treatment residual.
  dr <- crossfit_dr_pseudo_outcome(sim, "y_sim", cfg, seed_offset = seed_offset)
  x <- dr$data
  e <- x$propensity
  mu <- e * x$m1_hat + (1 - e) * x$m0_hat
  yres <- x$y_sim - mu
  ares <- x$treatment - e
  w <- normalize_weights(x$analysis_weight)
  s <- x$shock_sales_decline
  boundary <- range(s)
  knots <- weighted_quantile(s, w, unlist(cfg$second_stage$spline_knot_probs))
  b <- make_spline_basis(s, knots, boundary, cfg$second_stage$spline_degree)
  bg <- make_spline_basis(shock_grid, knots, boundary, cfg$second_stage$spline_degree)
  xreg <- b * ares
  lambda <- cfg$second_stage$lambda_spline_grid[[4]] %||% 0.1
  pen <- lambda * second_difference_penalty(ncol(b))
  normal <- crossprod(xreg, w * xreg) + pen + diag(1e-10, ncol(b))
  beta <- solve(normal, crossprod(xreg, w * yres))
  est <- as.numeric(bg %*% beta)

  residual <- yres - as.numeric(xreg %*% beta)
  ainv <- solve(normal)
  h <- (xreg %*% ainv %*% t(bg)) * w
  clusters <- as.character(x$economy)
  ug <- unique(clusters)
  scores <- matrix(0, length(ug), length(shock_grid))
  for (j in seq_along(ug)) {
    idx <- which(clusters == ug[j])
    scores[j, ] <- colSums(h[idx, , drop = FALSE] * residual[idx])
  }
  se <- sqrt(length(ug) / (length(ug) - 1) * colSums(scores^2))
  z <- stats::qnorm(1 - (1 - cfg$inference$pointwise_level) / 2)
  data.table::data.table(shock = shock_grid, estimate = est,
                         lower = est - z * se, upper = est + z * se)
}

fit_grf_curve <- function(sim, cfg, shock_grid, seed_offset) {
  cov <- available_covariates(sim, include_shock = TRUE)
  prep <- fit_preprocessor(sim, cov$continuous, cov$categorical, "analysis_weight")
  proc <- apply_preprocessor(sim, prep)
  xmat <- stats::model.matrix(~ . - 1, data = proc)
  forest <- grf::causal_forest(
    X = xmat, Y = sim$y_sim, W = sim$treatment,
    sample.weights = sim$analysis_weight,
    clusters = as.factor(sim$economy),
    num.trees = cfg$validation$grf_num_trees,
    seed = cfg$seed + seed_offset
  )

  rows <- lapply(shock_grid, function(s) {
    q <- data.table::copy(sim)
    q[, shock_sales_decline := s]
    qp <- apply_preprocessor(q, prep)
    qm <- stats::model.matrix(~ . - 1, data = qp)
    missing <- setdiff(colnames(xmat), colnames(qm))
    if (length(missing)) qm <- cbind(qm, matrix(0, nrow(qm), length(missing), dimnames = list(NULL, missing)))
    qm <- qm[, colnames(xmat), drop = FALSE]
    pred <- predict(forest, newdata = qm, estimate.variance = TRUE)
    ww <- normalize_weights(q$analysis_weight)
    estimate <- weighted_mean(pred$predictions, ww)
    # GRF reports unit-level variance estimates; aggregate them with a weighted
    # diagonal approximation for the standardized curve.
    se <- sqrt(sum((ww / sum(ww))^2 * pmax(pred$variance.estimates, 0), na.rm = TRUE))
    data.table::data.table(shock = s, estimate = estimate,
                           lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)
  })
  data.table::rbindlist(rows)
}

causalpfn_model_matrices <- function(sim, shock_grid, cfg, seed_offset = 0L) {
  cov <- available_covariates(sim, include_shock = TRUE)
  prep <- fit_preprocessor(sim, cov$continuous, cov$categorical, "analysis_weight")
  proc <- apply_preprocessor(sim, prep)
  xtrain <- stats::model.matrix(~ . - 1, data = proc)
  train <- data.table::as.data.table(xtrain)
  train[, `:=`(treatment = sim$treatment, outcome = sim$y_sim)]

  # Exact standardization would duplicate every firm at every grid point. To keep
  # query size bounded, represent the weighted target distribution with a fixed-size
  # Monte Carlo draw per grid; the draw count is configurable in analysis.yml.
  draws <- as.integer(cfg$validation$causalpfn_query_draws_per_grid)
  prob <- sim$analysis_weight / sum(sim$analysis_weight)
  set.seed(cfg$seed + 800000L + seed_offset)
  sampled_index <- sample(seq_len(nrow(sim)), draws, replace = TRUE, prob = prob)

  queries <- lapply(seq_along(shock_grid), function(j) {
    q <- data.table::copy(sim[sampled_index])
    q[, shock_sales_decline := shock_grid[j]]
    qp <- apply_preprocessor(q, prep)
    qm <- stats::model.matrix(~ . - 1, data = qp)
    missing <- setdiff(colnames(xtrain), colnames(qm))
    if (length(missing)) qm <- cbind(qm, matrix(0, nrow(qm), length(missing), dimnames = list(NULL, missing)))
    qm <- qm[, colnames(xtrain), drop = FALSE]
    z <- data.table::as.data.table(qm)
    z[, `:=`(grid_id = j, analysis_weight = 1)]
    z
  })
  list(train = train, query = data.table::rbindlist(queries))
}

fit_causalpfn_curve <- function(sim, cfg, shock_grid, seed_offset, root = project_root()) {
  mm <- causalpfn_model_matrices(sim, shock_grid, cfg, seed_offset)
  temp <- file.path(root, "data", "derived", "causalpfn_temp")
  dir.create(temp, recursive = TRUE, showWarnings = FALSE)
  tag <- paste0("seed_", cfg$seed + seed_offset)
  train_path <- file.path(temp, paste0(tag, "_train.csv"))
  query_path <- file.path(temp, paste0(tag, "_query.csv"))
  out_path <- file.path(temp, paste0(tag, "_out.csv"))
  summary_path <- file.path(temp, paste0(tag, "_summary.csv"))
  data.table::fwrite(mm$train, train_path)
  data.table::fwrite(mm$query, query_path)

  py <- Sys.which("python")
  if (!nzchar(py)) py <- Sys.which("python3")
  if (!nzchar(py)) stop("Python interpreter not found for CausalPFN")
  bridge <- file.path(root, "python", "causalpfn_bridge.py")
  status <- system2(py, c(
    bridge, "--train", train_path, "--query", query_path,
    "--output", out_path, "--summary", summary_path,
    "--with-ci", "--ci-samples", cfg$validation$causalpfn_ci_samples,
    "--seed", cfg$seed + seed_offset
  ))
  if (status != 0L || !file.exists(summary_path)) stop("CausalPFN bridge failed")
  out <- data.table::fread(summary_path)
  out[, shock := shock_grid[grid_id]]
  out[, .(shock, estimate, lower, upper)]
}

run_validation_repetition <- function(primary_fit, mechanism, repetition, cfg,
                                      methods = c("grf", "causalpfn", "model_averaged_dml", "conditional_dr")) {
  seed <- cfg$seed + 100000L + match(mechanism, c("smooth_transition", "sharp_jump", "linear_change")) * 10000L + repetition
  sim <- generate_semisynthetic_data(primary_fit, mechanism, cfg, seed)
  pars <- attr(sim, "dgp_parameters")
  grid <- seq(cfg$sample$primary_curve_min, cfg$sample$primary_curve_max,
              by = cfg$sample$shock_grid_step)
  true_curve <- data.table::data.table(
    shock = grid,
    true_effect = true_effect_function(grid, mechanism, pars,
                                       c(cfg$sample$primary_curve_min, cfg$sample$primary_curve_max))
  )

  estimators <- list(
    grf = function() fit_grf_curve(sim, cfg, grid, repetition),
    causalpfn = function() fit_causalpfn_curve(sim, cfg, grid, repetition),
    model_averaged_dml = function() fit_model_averaged_dml_curve(sim, cfg, grid, repetition),
    conditional_dr = function() fit_primary_validation_method(sim, cfg, grid, repetition)
  )
  rows <- lapply(methods, function(method) {
    curve <- estimators[[method]]()
    met <- curve_metrics(curve, true_curve, mechanism, pars)
    met[, `:=`(mechanism = mechanism, method = method, repetition = repetition)]
    met
  })
  data.table::rbindlist(rows)
}

run_component_repetition <- function(primary_fit, repetition, cfg) {
  mechanism <- "smooth_transition"
  seed <- cfg$seed + 300000L + repetition
  sim <- generate_semisynthetic_data(primary_fit, mechanism, cfg, seed)
  pars <- attr(sim, "dgp_parameters")
  grid <- seq(cfg$sample$primary_curve_min, cfg$sample$primary_curve_max,
              by = cfg$sample$shock_grid_step)
  true_curve <- data.table::data.table(
    shock = grid,
    true_effect = true_effect_function(grid, mechanism, pars,
                                       c(cfg$sample$primary_curve_min, cfg$sample$primary_curve_max))
  )
  curves <- list(
    full_specification = fit_primary_validation_method(sim, cfg, grid, 4000L + repetition),
    without_covariate_standardization = fit_unstandardized_spline(sim, cfg, grid, 5000L + repetition),
    without_dr_orthogonalization = fit_no_orthogonalization(sim, cfg, grid, 6000L + repetition)
  )
  data.table::rbindlist(lapply(names(curves), function(method) {
    met <- curve_metrics(curves[[method]], true_curve, mechanism, pars)
    met[, `:=`(mechanism = mechanism, method = method, repetition = repetition)]
    met
  }))
}

summarize_validation <- function(raw) {
  raw[, .(
    standardized_rmse = mean(standardized_rmse, na.rm = TRUE),
    change_location_error = if (all(is.na(change_location_error))) NA_real_ else mean(change_location_error, na.rm = TRUE),
    pointwise_coverage = mean(pointwise_coverage, na.rm = TRUE)
  ), by = .(mechanism, method)]
}
