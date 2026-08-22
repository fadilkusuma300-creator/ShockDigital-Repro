# Robustness, placebo, and negative-control analyses.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "03_preprocess.R"))
source(file.path(project_root(), "R", "05_dr_estimator.R"))
source(file.path(project_root(), "R", "06_continuous_curve.R"))

clustered_lm_result <- function(fit, term, cluster) {
  vc <- sandwich::vcovCL(fit, cluster = cluster, type = "HC1")
  if (!term %in% names(stats::coef(fit))) return(NULL)
  est <- stats::coef(fit)[term]
  se <- sqrt(vc[term, term])
  ci <- ci_from_estimate(est, se)
  data.table::data.table(estimate = est, se = se, lower = ci[1], upper = ci[2])
}

fixed_effects_severe <- function(dt, outcome = "y_sales") {
  x <- prepare_raw_analysis_columns(dt)
  x <- x[severe_shock == 1L & is.finite(get(outcome))]
  cov <- second_stage_covariates(x)

  # Economy, industry, and survey round enter explicitly as fixed effects here,
  # using the same pre-treatment adjustment set as the primary analysis.
  fixed_cat <- unique(c(cov$categorical, "economy", "industry", "survey_round"))
  fixed_cat <- intersect(fixed_cat, names(x))
  fixed_cont <- intersect(c("shock_sales_decline", cov$continuous), names(x))
  prep <- fit_preprocessor(x, fixed_cont, fixed_cat, "analysis_weight")
  proc <- apply_preprocessor(x, prep)
  dat <- data.frame(
    y = safe_numeric(x[[outcome]]),
    treatment = x$treatment,
    proc,
    check.names = FALSE
  )
  fit <- stats::lm(y ~ treatment + ., data = dat, weights = x$analysis_weight)
  out <- clustered_lm_result(fit, "treatment", x$economy)
  if (is.null(out)) stop("Treatment coefficient absent from fixed-effects model")
  out[, method := "fixed_effects"]
  out
}

strict_support_effect <- function(dr_fit, cfg) {
  lo <- cfg$sample$strict_propensity_lower
  hi <- cfg$sample$strict_propensity_upper
  x <- dr_fit$data[propensity_raw >= lo & propensity_raw <= hi & severe_shock == 1L]
  z <- cluster_mean_se(x$dr_pseudo_outcome, x$analysis_weight, x$economy)
  ci <- ci_from_estimate(z$estimate, z$se)
  data.table::data.table(
    method = "strict_support_0.10_0.90", estimate = z$estimate, se = z$se,
    lower = ci[1], upper = ci[2], n = z$n
  )
}

fit_attrition_probabilities <- function(all_windows, outcome, cfg) {
  # The paper specifies inverse-probability attrition correction but does not
  # disclose the exact attrition learner. We use a weighted logistic model with
  # the same pre-treatment W set plus treatment. This choice is isolated here so
  # an author-supplied specification can be substituted without other changes.
  x <- prepare_raw_analysis_columns(all_windows)
  x[, observed_outcome := as.integer(is.finite(get(outcome)))]
  cov <- available_covariates(x, include_shock = TRUE)
  prep <- fit_preprocessor(x, cov$continuous, cov$categorical, "analysis_weight")
  proc <- apply_preprocessor(x, prep)
  dat <- data.frame(observed_outcome = x$observed_outcome,
                    treatment = x$treatment, proc, check.names = FALSE)
  fit <- stats::glm(
    observed_outcome ~ ., data = dat, family = stats::binomial(),
    weights = x$analysis_weight
  )
  p <- as.numeric(stats::predict(fit, type = "response"))
  clip(p, 0.02, 0.98)
}

attrition_corrected_effect <- function(windows, dr_fit, outcome, cfg) {
  pobs <- fit_attrition_probabilities(windows, outcome, cfg)
  keys <- windows[, .(economy, firm_id, p_observed = pobs)]
  x <- merge(dr_fit$data, keys, by = c("economy", "firm_id"), all.x = TRUE)
  x <- x[severe_shock == 1L & is.finite(p_observed)]
  x[, attrition_weight := analysis_weight / p_observed]
  z <- cluster_mean_se(x$dr_pseudo_outcome, x$attrition_weight, x$economy)
  ci <- ci_from_estimate(z$estimate, z$se)
  data.table::data.table(
    method = "inverse_probability_attrition", estimate = z$estimate, se = z$se,
    lower = ci[1], upper = ci[2], n = z$n
  )
}

drop_one_economy <- function(windows, outcome, cfg) {
  economies <- sort(unique(as.character(windows$economy)))
  data.table::rbindlist(lapply(seq_along(economies), function(j) {
    econ <- economies[j]
    d <- windows[economy != econ]
    fit <- crossfit_dr_pseudo_outcome(d, outcome, cfg, seed_offset = 100L + j)
    z <- weighted_dr_effect(fit$data, quote(severe_shock == 1L))
    ci <- ci_from_estimate(z$estimate, z$se)
    data.table::data.table(
      omitted_economy = econ, estimate = z$estimate, se = z$se,
      lower = ci[1], upper = ci[2], n = z$n
    )
  }))
}

placebo_future_treatment_past_sales <- function(windows) {
  # The placebo outcome precedes t1 treatment. Because shock severity is derived
  # from this same t0 sales outcome, it is intentionally excluded from the
  # adjustment set to avoid conditioning directly on the placebo outcome.
  x <- prepare_raw_analysis_columns(windows)
  x <- x[is.finite(t0_sales_change)]
  cov <- available_covariates(x, include_shock = FALSE)
  fixed_cat <- unique(c(cov$categorical, "economy", "industry", "survey_round"))
  fixed_cat <- intersect(fixed_cat, names(x))
  prep <- fit_preprocessor(x, cov$continuous, fixed_cat, "analysis_weight")
  proc <- apply_preprocessor(x, prep)
  dat <- data.frame(y = x$t0_sales_change, treatment = x$treatment, proc, check.names = FALSE)
  fit <- stats::lm(y ~ treatment + ., data = dat, weights = x$analysis_weight)
  ans <- clustered_lm_result(fit, "treatment", x$economy)
  ans[, method := "temporal_placebo"]
  ans
}

negative_control_export <- function(windows) {
  # Pre-pandemic export status is temporally fixed before the digital adjustment.
  # Exporter itself is removed from the adjustment set because it is the outcome.
  x <- prepare_raw_analysis_columns(windows)
  x <- x[!is.na(exporter)]
  cov <- available_covariates(x, include_shock = TRUE)
  cov$categorical <- setdiff(cov$categorical, "exporter")
  prep <- fit_preprocessor(x, cov$continuous, cov$categorical, "analysis_weight")
  proc <- apply_preprocessor(x, prep)
  dat <- data.frame(y = safe_numeric(x$exporter), treatment = x$treatment, proc,
                    check.names = FALSE)
  fit <- stats::lm(y ~ treatment + ., data = dat, weights = x$analysis_weight)
  ans <- clustered_lm_result(fit, "treatment", x$economy)
  ans[, `:=`(estimate = 100 * estimate, se = 100 * se,
             lower = 100 * lower, upper = 100 * upper,
             method = "negative_control_export")]
  ans
}

run_sales_robustness <- function(windows, primary_fit, cfg) {
  rows <- list()

  main <- weighted_dr_effect(primary_fit$data, quote(severe_shock == 1L))
  ci <- ci_from_estimate(main$estimate, main$se)
  rows[["main"]] <- data.table::data.table(
    method = "conditional_dr", estimate = main$estimate, se = main$se,
    lower = ci[1], upper = ci[2], n = main$n
  )

  rows[["fixed"]] <- fixed_effects_severe(windows, "y_sales")

  ipw <- ipw_effect(primary_fit$data[severe_shock == 1L], "y_sales")
  ci <- ci_from_estimate(ipw$estimate, ipw$se)
  rows[["ipw"]] <- data.table::data.table(
    method = "ipw", estimate = ipw$estimate, se = ipw$se,
    lower = ci[1], upper = ci[2], n = ipw$n
  )

  rows[["strict"]] <- strict_support_effect(primary_fit, cfg)
  rows[["attrition"]] <- attrition_corrected_effect(windows, primary_fit, "y_sales", cfg)

  alt <- crossfit_dr_pseudo_outcome(
    windows, "y_exit_severe_sales_decline", cfg, seed_offset = 301L
  )
  z <- weighted_dr_effect(alt$data, quote(severe_shock == 1L))
  ci <- ci_from_estimate(z$estimate, z$se)
  rows[["alt_sales"]] <- data.table::data.table(
    method = "alternative_sales_recovery", estimate = 100 * z$estimate, se = 100 * z$se,
    lower = 100 * ci[1], upper = 100 * ci[2], n = z$n
  )

  no50 <- windows[shock_sales_decline != cfg$sample$severe_shock_cutoff]
  fit_no50 <- crossfit_dr_pseudo_outcome(no50, "y_sales", cfg, seed_offset = 302L)
  z <- weighted_dr_effect(fit_no50$data, quote(severe_shock == 1L))
  ci <- ci_from_estimate(z$estimate, z$se)
  rows[["exclude50"]] <- data.table::data.table(
    method = "exclude_exact_50", estimate = z$estimate, se = z$se,
    lower = ci[1], upper = ci[2], n = z$n
  )

  strict_input <- crossfit_dr_pseudo_outcome(windows, "y_input_increase", cfg, seed_offset = 303L)
  z <- weighted_dr_effect(strict_input$data, quote(severe_shock == 1L))
  ci <- ci_from_estimate(z$estimate, z$se)
  rows[["strict_input"]] <- data.table::data.table(
    method = "strict_input_recovery", estimate = 100 * z$estimate, se = 100 * z$se,
    lower = 100 * ci[1], upper = 100 * ci[2], n = z$n
  )

  rows[["placebo"]] <- placebo_future_treatment_past_sales(windows)
  rows[["negative"]] <- negative_control_export(windows)

  list(
    summary = data.table::rbindlist(rows, fill = TRUE),
    drop_one_economy = drop_one_economy(windows, "y_sales", cfg)
  )
}
