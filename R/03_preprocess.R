# Fold-specific preprocessing for nuisance models and balance diagnostics.

source(file.path(project_root(), "R", "00_utils.R"))

analysis_covariates <- function() {
  # S is kept separate so the same Z list can be used for standardization.
  list(
    continuous = c(
      "baseline_employment", "firm_age", "manager_experience",
      "relative_interview_month"
    ),
    categorical = c(
      "foreign_owned", "exporter", "website", "credit_constraints", "innovation",
      "t0_demand_state", "t0_cash_flow_state", "t0_input_state",
      "t0_operating_status", "industry", "economy", "survey_round"
    )
  )
}

available_covariates <- function(dt, include_shock = FALSE) {
  spec <- analysis_covariates()
  cont <- intersect(spec$continuous, names(dt))
  cat <- intersect(spec$categorical, names(dt))
  if (include_shock) cont <- c("shock_sales_decline", cont)
  list(continuous = unique(cont), categorical = unique(cat))
}

fit_preprocessor <- function(train, continuous, categorical, weight_col = "analysis_weight") {
  w <- train[[weight_col]]
  prep <- list(continuous = continuous, categorical = categorical,
               medians = list(), modes = list(), levels = list())

  for (nm in continuous) {
    x <- safe_numeric(train[[nm]])
    med <- weighted_quantile(x, w, 0.5)[1]
    if (!is.finite(med)) med <- 0
    prep$medians[[nm]] <- med
  }

  for (nm in categorical) {
    x <- safe_character(train[[nm]])
    x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
    mode_value <- weighted_mode(x, w)
    if (is.na(mode_value) || !nzchar(mode_value)) mode_value <- "__MISSING__"
    observed <- sort(unique(x[!is.na(x)]))
    prep$modes[[nm]] <- as.character(mode_value)
    prep$levels[[nm]] <- unique(c(observed, as.character(mode_value), "__OTHER__"))
  }
  class(prep) <- "survey_preprocessor"
  prep
}

apply_preprocessor <- function(dt, prep) {
  out <- data.table::data.table(row_id = seq_len(nrow(dt)))

  for (nm in prep$continuous) {
    x <- safe_numeric(dt[[nm]])
    miss <- as.integer(!is.finite(x))
    x[!is.finite(x)] <- prep$medians[[nm]]
    out[[nm]] <- x
    out[[paste0(nm, "__missing")]] <- miss
  }

  for (nm in prep$categorical) {
    x <- safe_character(dt[[nm]])
    x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
    miss <- is.na(x)
    x[miss] <- prep$modes[[nm]]
    known <- prep$levels[[nm]]
    x[!x %in% known] <- "__OTHER__"
    out[[nm]] <- factor(x, levels = known)
    out[[paste0(nm, "__missing")]] <- as.integer(miss)
  }
  out[, row_id := NULL]
  out
}

make_model_matrix_pair <- function(train, test, continuous, categorical,
                                   weight_col = "analysis_weight") {
  prep <- fit_preprocessor(train, continuous, categorical, weight_col)
  tr <- apply_preprocessor(train, prep)
  te <- apply_preprocessor(test, prep)

  # Build both matrices with the same factor levels. model.matrix may still drop
  # absent columns, so explicit alignment is performed after construction.
  f <- stats::as.formula("~ . - 1")
  xtr <- stats::model.matrix(f, data = tr)
  xte <- stats::model.matrix(f, data = te)
  all_cols <- union(colnames(xtr), colnames(xte))

  align <- function(x, cols) {
    missing <- setdiff(cols, colnames(x))
    if (length(missing)) {
      x <- cbind(x, matrix(0, nrow(x), length(missing),
                           dimnames = list(NULL, missing)))
    }
    x[, cols, drop = FALSE]
  }
  xtr <- align(xtr, all_cols)
  xte <- align(xte, all_cols)

  # Remove constant columns based only on the training fold.
  variable <- apply(xtr, 2, function(z) {
    z <- z[is.finite(z)]
    length(unique(z)) > 1L
  })
  if (!any(variable)) stop("No varying nuisance-model covariates after preprocessing")
  xtr <- xtr[, variable, drop = FALSE]
  xte <- xte[, colnames(xtr), drop = FALSE]

  list(train = xtr, test = xte, preprocessor = prep)
}

prepare_raw_analysis_columns <- function(dt) {
  x <- data.table::copy(dt)

  # Recode common yes/no questionnaire fields to 0/1 when they retain public
  # World Bank coding (1=yes, 2=no). Existing 0/1 indicators are left unchanged.
  for (nm in intersect(c("website", "innovation"), names(x))) {
    z <- safe_numeric(x[[nm]])
    vals <- unique(z[is.finite(z)])
    if (length(vals) && all(vals %in% c(1, 2))) {
      x[, (nm) := data.table::fcase(z == 1, 1L, z == 2, 0L, default = NA_integer_)]
    }
  }

  # Credit constraints are harmonized upstream to a binary indicator using the
  # World Bank fin23/fin24/fin25 classification (or a validated direct source).
  if ("credit_constraints" %in% names(x)) {
    zc <- safe_numeric(x$credit_constraints)
    bad <- is.finite(zc) & !zc %in% c(0, 1)
    if (any(bad)) stop("credit_constraints must be binary 0/1 after harmonization")
    x[, credit_constraints := as.integer(zc)]
  }
  x[, shock_sales_decline := pmax(0, pmin(100, safe_numeric(shock_sales_decline)))]
  x
}
