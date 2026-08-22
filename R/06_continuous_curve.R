# Penalized continuous-effect second stage.
# The model combines a cubic B-spline in shock severity, fixed adjustment terms,
# and ridge-shrunk economy/industry intercepts. For a Gaussian second stage,
# ridge penalties on group intercepts are the penalized least-squares form of
# hierarchical random-intercept shrinkage.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "03_preprocess.R"))

second_stage_covariates <- function(dt) {
  spec <- analysis_covariates()
  list(
    continuous = intersect(spec$continuous, names(dt)),
    categorical = setdiff(intersect(spec$categorical, names(dt)), c("economy", "industry"))
  )
}

make_spline_basis <- function(s, knots, boundary, degree = 3L) {
  # Repeated/degenerate knots can occur in highly concentrated shock data.
  knots <- sort(unique(knots[is.finite(knots) & knots > boundary[1] & knots < boundary[2]]))
  splines::bs(
    s, knots = knots, degree = degree, Boundary.knots = boundary,
    intercept = FALSE, warn.outside = FALSE
  )
}

fit_second_stage_schema <- function(train, cfg) {
  cov <- second_stage_covariates(train)
  prep <- fit_preprocessor(train, cov$continuous, cov$categorical, "analysis_weight")
  proc <- apply_preprocessor(train, prep)
  fixed <- stats::model.matrix(~ ., data = proc)

  # Keep the intercept and remove redundant zero-variance adjustment columns.
  keep <- vapply(seq_len(ncol(fixed)), function(j) {
    if (colnames(fixed)[j] == "(Intercept)") return(TRUE)
    length(unique(fixed[, j])) > 1L
  }, logical(1))
  fixed <- fixed[, keep, drop = FALSE]

  s <- train$shock_sales_decline
  w <- train$analysis_weight
  boundary <- range(s[is.finite(s)], na.rm = TRUE)
  if (diff(boundary) <= 0) stop("Shock severity has no variation")
  knots <- weighted_quantile(s, w, unlist(cfg$second_stage$spline_knot_probs))
  knots <- sort(unique(knots[is.finite(knots) & knots > boundary[1] & knots < boundary[2]]))
  spline <- make_spline_basis(s, knots, boundary, cfg$second_stage$spline_degree)
  colnames(spline) <- paste0("spline_", seq_len(ncol(spline)))

  economy_levels <- sort(unique(as.character(train$economy)))
  industry_levels <- sort(unique(as.character(train$industry)))

  list(
    preprocessor = prep,
    fixed_columns = colnames(fixed),
    knots = knots,
    boundary = boundary,
    degree = cfg$second_stage$spline_degree,
    economy_levels = economy_levels,
    industry_levels = industry_levels,
    n_fixed = ncol(fixed),
    n_spline = ncol(spline),
    n_economy = length(economy_levels),
    n_industry = length(industry_levels)
  )
}

align_fixed_matrix <- function(newdata, schema) {
  proc <- apply_preprocessor(newdata, schema$preprocessor)
  x <- stats::model.matrix(~ ., data = proc)
  missing <- setdiff(schema$fixed_columns, colnames(x))
  if (length(missing)) {
    x <- cbind(x, matrix(0, nrow(x), length(missing), dimnames = list(NULL, missing)))
  }
  x <- x[, schema$fixed_columns, drop = FALSE]
  x
}

random_dummy_matrix <- function(values, levels, prefix) {
  n <- length(values)
  m <- matrix(0, n, length(levels))
  colnames(m) <- paste0(prefix, levels)
  idx <- match(as.character(values), levels)
  ok <- which(!is.na(idx))
  if (length(ok)) m[cbind(ok, idx[ok])] <- 1
  m
}

build_second_stage_design <- function(newdata, schema, shock_override = NULL) {
  fixed <- align_fixed_matrix(newdata, schema)
  s <- if (is.null(shock_override)) newdata$shock_sales_decline else rep(shock_override, nrow(newdata))
  spline <- make_spline_basis(s, schema$knots, schema$boundary, schema$degree)
  colnames(spline) <- paste0("spline_", seq_len(ncol(spline)))
  econ <- random_dummy_matrix(newdata$economy, schema$economy_levels, "economy::")
  ind <- random_dummy_matrix(newdata$industry, schema$industry_levels, "industry::")
  cbind(fixed, spline, econ, ind)
}

second_difference_penalty <- function(k) {
  if (k < 3L) return(diag(k))
  d <- diff(diag(k), differences = 2L)
  crossprod(d)
}

second_stage_penalty <- function(schema, lambda_spline, lambda_economy, lambda_industry) {
  p <- schema$n_fixed + schema$n_spline + schema$n_economy + schema$n_industry
  pen <- matrix(0, p, p)
  i_s <- schema$n_fixed + seq_len(schema$n_spline)
  i_e <- schema$n_fixed + schema$n_spline + seq_len(schema$n_economy)
  i_i <- schema$n_fixed + schema$n_spline + schema$n_economy + seq_len(schema$n_industry)
  if (length(i_s)) pen[i_s, i_s] <- lambda_spline * second_difference_penalty(schema$n_spline)
  if (length(i_e)) pen[i_e, i_e] <- lambda_economy * diag(length(i_e))
  if (length(i_i)) pen[i_i, i_i] <- lambda_industry * diag(length(i_i))
  pen
}

fit_penalized_second_stage <- function(train, cfg, lambda_spline, lambda_economy,
                                       lambda_industry, schema = NULL) {
  if (is.null(schema)) schema <- fit_second_stage_schema(train, cfg)
  x <- build_second_stage_design(train, schema)
  y <- train$dr_pseudo_outcome
  w <- normalize_weights(train$analysis_weight)
  pen <- second_stage_penalty(schema, lambda_spline, lambda_economy, lambda_industry)

  # Weighted penalized normal equations. A tiny numerical ridge is added only
  # for stable inversion; it is many orders below the tuned penalties.
  xtw <- t(x) * rep(w, each = ncol(x))
  normal <- xtw %*% x + pen + diag(1e-10, ncol(x))
  rhs <- xtw %*% y
  beta <- as.numeric(solve(normal, rhs))
  fitted <- as.numeric(x %*% beta)

  structure(list(
    beta = beta, schema = schema, normal = normal,
    normal_inv = solve(normal), design = x, weights = w,
    fitted = fitted, residuals = y - fitted,
    lambda_spline = lambda_spline,
    lambda_economy = lambda_economy,
    lambda_industry = lambda_industry
  ), class = "penalized_effect_curve")
}

predict_penalized_second_stage <- function(model, newdata, shock_override = NULL) {
  x <- build_second_stage_design(newdata, model$schema, shock_override)
  as.numeric(x %*% model$beta)
}

second_stage_cv_score <- function(train, cfg, lambdas, folds) {
  loss <- numeric(max(folds))
  for (f in seq_len(max(folds))) {
    tr <- train[folds != f]
    va <- train[folds == f]
    if (!nrow(tr) || !nrow(va)) return(Inf)
    schema <- tryCatch(fit_second_stage_schema(tr, cfg), error = function(e) NULL)
    if (is.null(schema)) return(Inf)
    fit <- tryCatch(
      fit_penalized_second_stage(
        tr, cfg, lambdas[1], lambdas[2], lambdas[3], schema
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(Inf)
    pred <- predict_penalized_second_stage(fit, va)
    loss[f] <- weighted_mean((va$dr_pseudo_outcome - pred)^2, va$analysis_weight)
  }
  weighted_mean(loss, rep(1, length(loss)))
}

tune_second_stage <- function(train, cfg, seed_offset = 0L) {
  folds <- make_group_folds(
    train$economy, cfg$second_stage$inner_folds,
    cfg$seed + 5000L + seed_offset
  )
  grid <- data.table::CJ(
    lambda_spline = unlist(cfg$second_stage$lambda_spline_grid),
    lambda_economy = unlist(cfg$second_stage$lambda_economy_grid),
    lambda_industry = unlist(cfg$second_stage$lambda_industry_grid),
    sorted = FALSE
  )
  grid[, weighted_mse := vapply(seq_len(.N), function(i) {
    second_stage_cv_score(
      train, cfg,
      c(lambda_spline[i], lambda_economy[i], lambda_industry[i]), folds
    )
  }, numeric(1))]
  data.table::setorder(grid, weighted_mse, lambda_spline, lambda_economy, lambda_industry)
  grid
}

fit_continuous_effect_curve <- function(dr_data, cfg, seed_offset = 0L) {
  x <- data.table::copy(dr_data)
  grid_cv <- tune_second_stage(x, cfg, seed_offset)
  best <- grid_cv[1]
  fit <- fit_penalized_second_stage(
    x, cfg, best$lambda_spline, best$lambda_economy, best$lambda_industry
  )
  fit$tuning_grid <- grid_cv
  fit$data <- x
  fit
}

standardization_rows <- function(model, target_data, shock_grid) {
  # c(s) is the survey-weighted average design vector after setting every target
  # establishment to the same shock value s while preserving its Z values.
  w <- normalize_weights(target_data$analysis_weight)
  out <- matrix(NA_real_, length(shock_grid), length(model$beta))
  for (j in seq_along(shock_grid)) {
    d <- build_second_stage_design(target_data, model$schema, shock_grid[j])
    out[j, ] <- colSums(d * w) / sum(w)
  }
  rownames(out) <- as.character(shock_grid)
  out
}

standardized_effect_curve <- function(model, target_data, shock_grid) {
  cmat <- standardization_rows(model, target_data, shock_grid)
  data.table::data.table(
    shock = shock_grid,
    estimate = as.numeric(cmat %*% model$beta)
  )
}
