# Cross-fitted nuisance models for propensity and potential outcomes.
# Outcome regression uses a convex validation-weighted ensemble of elastic net,
# gradient boosting, and ridge models as described in the study.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "03_preprocess.R"))

fit_propensity_rf <- function(x, a, w, cfg) {
  p <- ncol(x)
  mtry <- max(1L, floor(cfg$nuisance_models$propensity$mtry_fraction * p))
  dat <- data.frame(a = factor(a, levels = c(0, 1)), x, check.names = FALSE)
  ranger::ranger(
    dependent.variable.name = "a",
    data = dat,
    probability = TRUE,
    num.trees = cfg$nuisance_models$propensity$num_trees,
    min.node.size = cfg$nuisance_models$propensity$min_node_size,
    mtry = mtry,
    case.weights = w,
    importance = "none",
    seed = cfg$seed,
    num.threads = cfg$runtime$parallel_cores
  )
}

predict_propensity_rf <- function(model, x) {
  p <- predict(model, data = data.frame(x, check.names = FALSE))$predictions
  if (is.matrix(p)) {
    if ("1" %in% colnames(p)) return(as.numeric(p[, "1"]))
    return(as.numeric(p[, ncol(p)]))
  }
  as.numeric(p)
}

weighted_train_validation_split <- function(dt, fraction, seed) {
  set.seed(seed)
  # Stratify by treatment and economy where possible to stabilize ensemble
  # weight selection inside each outer training fold.
  strata <- interaction(dt$economy, dt$treatment, drop = TRUE)
  val <- logical(nrow(dt))
  for (lev in levels(strata)) {
    idx <- which(strata == lev)
    nval <- max(1L, floor(length(idx) * fraction))
    if (length(idx) > 1L) nval <- min(nval, length(idx) - 1L)
    if (length(idx) == 1L) nval <- 0L
    if (nval > 0L) val[sample(idx, nval)] <- TRUE
  }
  if (!any(val)) val[sample(seq_len(nrow(dt)), max(1L, floor(nrow(dt) * fraction)))] <- TRUE
  list(train = which(!val), validation = which(val))
}

fit_glmnet_regression <- function(x, y, w, alpha, nlambda, seed) {
  set.seed(seed)
  # Internal lambda selection occurs exclusively inside the current training
  # fold. Fold IDs are supplied explicitly for deterministic weighted CV.
  nfolds <- min(5L, max(3L, floor(length(y) / 25L)))
  foldid <- rep(seq_len(nfolds), length.out = length(y))
  foldid <- sample(foldid, length(foldid), replace = FALSE)
  glmnet::cv.glmnet(
    x = x, y = y, weights = w, family = "gaussian", alpha = alpha,
    nlambda = nlambda, foldid = foldid, standardize = TRUE,
    intercept = TRUE, type.measure = "mse"
  )
}

predict_glmnet_regression <- function(model, x) {
  as.numeric(stats::predict(model, newx = x, s = "lambda.min"))
}

fit_xgboost_regression <- function(x_train, y_train, w_train,
                                   x_val = NULL, y_val = NULL, w_val = NULL,
                                   cfg, seed) {
  pars <- list(
    objective = "reg:squarederror",
    eta = cfg$nuisance_models$xgboost$eta,
    max_depth = cfg$nuisance_models$xgboost$max_depth,
    min_child_weight = cfg$nuisance_models$xgboost$min_child_weight,
    subsample = cfg$nuisance_models$xgboost$subsample,
    colsample_bytree = cfg$nuisance_models$xgboost$colsample_bytree,
    eval_metric = "rmse",
    seed = seed,
    nthread = cfg$runtime$parallel_cores
  )
  dtrain <- xgboost::xgb.DMatrix(x_train, label = y_train, weight = w_train)
  watch <- list(train = dtrain)
  early <- NULL
  if (!is.null(x_val) && nrow(x_val) > 0L) {
    dval <- xgboost::xgb.DMatrix(x_val, label = y_val, weight = w_val)
    watch$validation <- dval
    early <- cfg$nuisance_models$xgboost$early_stopping_rounds
  }
  xgboost::xgb.train(
    params = pars,
    data = dtrain,
    nrounds = cfg$nuisance_models$xgboost$nrounds,
    watchlist = watch,
    early_stopping_rounds = early,
    verbose = 0
  )
}

predict_xgboost_regression <- function(model, x) {
  as.numeric(stats::predict(model, xgboost::xgb.DMatrix(x)))
}

solve_convex_ensemble <- function(pred_matrix, y, w) {
  ok <- is.finite(y) & is.finite(w) & w > 0 & apply(pred_matrix, 1, function(r) all(is.finite(r)))
  p <- pred_matrix[ok, , drop = FALSE]
  y <- y[ok]
  w <- w[ok]
  k <- ncol(p)
  if (nrow(p) < k + 2L) return(rep(1 / k, k))

  sw <- sqrt(w / mean(w))
  pw <- p * sw
  yw <- y * sw
  dmat <- 2 * crossprod(pw) + diag(1e-8, k)
  dvec <- 2 * crossprod(pw, yw)[, 1]
  # solve.QP uses A' b >= b0. First constraint is equality sum(weights)=1;
  # remaining constraints enforce non-negative ensemble weights.
  amat <- cbind(rep(1, k), diag(k))
  bvec <- c(1, rep(0, k))
  fit <- tryCatch(
    quadprog::solve.QP(dmat, dvec, amat, bvec, meq = 1L),
    error = function(e) NULL
  )
  if (is.null(fit)) return(rep(1 / k, k))
  alpha <- pmax(0, fit$solution)
  if (sum(alpha) <= 0) return(rep(1 / k, k))
  alpha / sum(alpha)
}

fit_outcome_ensemble <- function(train_dt, outcome, covariates, cfg, seed) {
  split <- weighted_train_validation_split(
    train_dt, cfg$cross_fitting$inner_validation_fraction, seed
  )
  tr <- train_dt[split$train]
  va <- train_dt[split$validation]

  mm <- make_model_matrix_pair(
    tr, va, covariates$continuous, covariates$categorical, "analysis_weight"
  )
  xtr <- mm$train
  xva <- mm$test
  ytr <- safe_numeric(tr[[outcome]])
  yva <- safe_numeric(va[[outcome]])
  wtr <- tr$analysis_weight
  wva <- va$analysis_weight

  en <- fit_glmnet_regression(
    xtr, ytr, wtr, cfg$nuisance_models$elastic_net$alpha,
    cfg$nuisance_models$elastic_net$nlambda, seed + 1L
  )
  ridge <- fit_glmnet_regression(
    xtr, ytr, wtr, cfg$nuisance_models$ridge$alpha,
    cfg$nuisance_models$ridge$nlambda, seed + 2L
  )
  boost <- fit_xgboost_regression(
    xtr, ytr, wtr, xva, yva, wva, cfg, seed + 3L
  )

  pval <- cbind(
    elastic_net = predict_glmnet_regression(en, xva),
    gradient_boosting = predict_xgboost_regression(boost, xva),
    ridge = predict_glmnet_regression(ridge, xva)
  )
  alpha <- solve_convex_ensemble(pval, yva, wva)
  names(alpha) <- colnames(pval)

  # Refit all learners on the complete outer-training arm. Preprocessing is also
  # refit so no outer test information enters imputation or factor handling.
  prep_full <- fit_preprocessor(
    train_dt, covariates$continuous, covariates$categorical, "analysis_weight"
  )
  full_proc <- apply_preprocessor(train_dt, prep_full)
  xfull <- stats::model.matrix(~ . - 1, data = full_proc)
  variable <- apply(xfull, 2, function(z) length(unique(z[is.finite(z)])) > 1L)
  xfull <- xfull[, variable, drop = FALSE]

  en_full <- fit_glmnet_regression(
    xfull, safe_numeric(train_dt[[outcome]]), train_dt$analysis_weight,
    cfg$nuisance_models$elastic_net$alpha, cfg$nuisance_models$elastic_net$nlambda,
    seed + 11L
  )
  ridge_full <- fit_glmnet_regression(
    xfull, safe_numeric(train_dt[[outcome]]), train_dt$analysis_weight,
    cfg$nuisance_models$ridge$alpha, cfg$nuisance_models$ridge$nlambda,
    seed + 12L
  )

  # Use the validation-selected boosting iteration when available.
  cfg_refit <- cfg
  best_iter <- boost$best_iteration %||% cfg$nuisance_models$xgboost$nrounds
  cfg_refit$nuisance_models$xgboost$nrounds <- as.integer(best_iter)
  boost_full <- fit_xgboost_regression(
    xfull, safe_numeric(train_dt[[outcome]]), train_dt$analysis_weight,
    cfg = cfg_refit, seed = seed + 13L
  )

  structure(list(
    preprocessor = prep_full,
    columns = colnames(xfull),
    elastic_net = en_full,
    gradient_boosting = boost_full,
    ridge = ridge_full,
    ensemble_weights = alpha
  ), class = "outcome_ensemble")
}

matrix_for_ensemble <- function(model, newdata) {
  proc <- apply_preprocessor(newdata, model$preprocessor)
  x <- stats::model.matrix(~ . - 1, data = proc)
  missing <- setdiff(model$columns, colnames(x))
  if (length(missing)) {
    x <- cbind(x, matrix(0, nrow(x), length(missing), dimnames = list(NULL, missing)))
  }
  extra <- setdiff(colnames(x), model$columns)
  if (length(extra)) x <- x[, setdiff(colnames(x), extra), drop = FALSE]
  x[, model$columns, drop = FALSE]
}

predict_outcome_ensemble <- function(model, newdata) {
  x <- matrix_for_ensemble(model, newdata)
  p <- cbind(
    elastic_net = predict_glmnet_regression(model$elastic_net, x),
    gradient_boosting = predict_xgboost_regression(model$gradient_boosting, x),
    ridge = predict_glmnet_regression(model$ridge, x)
  )
  as.numeric(p %*% model$ensemble_weights[colnames(p)])
}
