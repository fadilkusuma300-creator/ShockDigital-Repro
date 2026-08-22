# Five-fold cross-fitted doubly robust estimation.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "03_preprocess.R"))
source(file.path(project_root(), "R", "04_nuisance_models.R"))

crossfit_dr_pseudo_outcome <- function(dt, outcome, cfg, support = c(0.02, 0.98), seed_offset = 0L) {
  x <- prepare_raw_analysis_columns(dt)
  needed <- c(outcome, "treatment", "analysis_weight", "economy", "shock_sales_decline")
  ok <- stats::complete.cases(x[, ..needed]) & is.finite(x[[outcome]])
  x <- x[ok]
  if (!nrow(x)) stop("No complete observations for outcome: ", outcome)
  if (length(unique(x$treatment)) < 2L) stop("Both treatment arms are required")

  z <- available_covariates(x, include_shock = TRUE)
  folds <- make_stratified_folds(
    x$economy, x$treatment,
    k = cfg$cross_fitting$outer_folds,
    seed = cfg$seed + seed_offset
  )

  e_hat <- m1_hat <- m0_hat <- rep(NA_real_, nrow(x))
  ensemble_weights <- list()

  for (fold in seq_len(cfg$cross_fitting$outer_folds)) {
    test_idx <- which(folds == fold)
    train_idx <- which(folds != fold)
    tr <- x[train_idx]
    te <- x[test_idx]

    # Propensity preprocessing is trained only on the outer training fold.
    mm <- make_model_matrix_pair(
      tr, te, z$continuous, z$categorical, "analysis_weight"
    )
    prop <- fit_propensity_rf(mm$train, tr$treatment, tr$analysis_weight, cfg)
    e_hat[test_idx] <- predict_propensity_rf(prop, mm$test)

    # Separate arm-specific outcome regressions m1(W) and m0(W).
    tr1 <- tr[treatment == 1L]
    tr0 <- tr[treatment == 0L]
    if (nrow(tr1) < 30L || nrow(tr0) < 30L) {
      stop("Too few observations in one treatment arm inside outer fold ", fold)
    }
    mod1 <- fit_outcome_ensemble(tr1, outcome, z, cfg, cfg$seed + 1000L + fold + seed_offset)
    mod0 <- fit_outcome_ensemble(tr0, outcome, z, cfg, cfg$seed + 2000L + fold + seed_offset)
    m1_hat[test_idx] <- predict_outcome_ensemble(mod1, te)
    m0_hat[test_idx] <- predict_outcome_ensemble(mod0, te)
    ensemble_weights[[fold]] <- data.table::data.table(
      fold = fold,
      arm = rep(c("treated", "control"), each = 3L),
      learner = rep(names(mod1$ensemble_weights), 2L),
      weight = c(mod1$ensemble_weights, mod0$ensemble_weights)
    )
  }

  e_raw <- e_hat
  in_support <- is.finite(e_raw) & e_raw >= support[1] & e_raw <= support[2]
  # Common-support observations are retained; observations outside the stated
  # interval are trimmed rather than moved onto the boundary. A clipped helper
  # is used only to keep the pre-trim arithmetic finite. For every retained row
  # it is numerically identical to the raw propensity score.
  e_used <- clip(e_raw, support[1], support[2])
  a <- x$treatment
  y <- safe_numeric(x[[outcome]])
  psi <- m1_hat - m0_hat + a * (y - m1_hat) / e_used -
    (1 - a) * (y - m0_hat) / (1 - e_used)

  x[, `:=`(
    outer_fold = folds,
    propensity_raw = e_raw,
    propensity = e_used,
    in_common_support = in_support,
    m1_hat = m1_hat,
    m0_hat = m0_hat,
    dr_pseudo_outcome = psi
  )]
  x <- x[in_common_support == TRUE]

  list(
    data = x,
    ensemble_weights = data.table::rbindlist(ensemble_weights),
    support = support,
    outcome = outcome
  )
}

weighted_dr_effect <- function(dr_data, subgroup = NULL) {
  x <- dr_data
  if (!is.null(subgroup)) x <- x[eval(subgroup)]
  cluster_mean_se(x$dr_pseudo_outcome, x$analysis_weight, x$economy)
}

estimate_ate_table <- function(dr_fit, cfg) {
  all_eff <- weighted_dr_effect(dr_fit$data)
  sev_eff <- weighted_dr_effect(dr_fit$data, quote(severe_shock == 1L))
  mk <- function(label, z) {
    ci <- ci_from_estimate(z$estimate, z$se, cfg$inference$pointwise_level)
    data.table::data.table(
      estimand = label, estimate = z$estimate, se = z$se,
      lower = ci[1], upper = ci[2], n = z$n, clusters = z$clusters
    )
  }
  data.table::rbindlist(list(mk("full_sample_ate", all_eff), mk("severe_shock_ate", sev_eff)))
}

ipw_effect <- function(dt, outcome, support = c(0.02, 0.98)) {
  x <- dt[is.finite(get(outcome)) & is.finite(propensity_raw)]
  x[, e_ipw := clip(propensity_raw, support[1], support[2])]
  x[, ipw_signal := treatment * get(outcome) / e_ipw -
        (1 - treatment) * get(outcome) / (1 - e_ipw)]
  cluster_mean_se(x$ipw_signal, x$analysis_weight, x$economy)
}
