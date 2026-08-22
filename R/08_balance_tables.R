# Covariate balance before and after inverse-propensity weighting.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "03_preprocess.R"))

binary_or_numeric_smd <- function(x, a, w) {
  if (is.factor(x) || is.character(x)) return(NA_real_)
  z <- safe_numeric(x)
  vals <- sort(unique(z[is.finite(z)]))
  if (length(vals) <= 2L && all(vals %in% c(0, 1))) {
    weighted_smd_binary(z, a, w)
  } else {
    weighted_smd_numeric(z, a, w)
  }
}

categorical_level_smds <- function(x, a, w, prefix) {
  z <- safe_character(x)
  z[is.na(z) | z %in% c("", "NA", "NaN", "NULL")] <- "__MISSING__"
  lev <- sort(unique(z))
  if (!length(lev)) return(data.table::data.table())
  data.table::rbindlist(lapply(lev, function(l) {
    b <- as.integer(z == l)
    data.table::data.table(
      covariate = paste0(prefix, "=", l),
      smd = weighted_smd_binary(b, a, w)
    )
  }))
}

balance_table <- function(dr_data) {
  x <- dr_data
  spec <- available_covariates(x, include_shock = TRUE)
  a <- x$treatment
  w0 <- x$analysis_weight
  e <- clip(x$propensity, 0.02, 0.98)
  # Stabilized survey-weighted ATE weights.
  wp <- w0 * (a / e + (1 - a) / (1 - e))

  rows <- list()
  for (nm in spec$continuous) {
    rows[[length(rows) + 1L]] <- data.table::data.table(
      covariate = nm,
      smd_before = binary_or_numeric_smd(x[[nm]], a, w0),
      smd_after = binary_or_numeric_smd(x[[nm]], a, wp)
    )
    miss <- as.integer(!is.finite(safe_numeric(x[[nm]])))
    if (any(miss == 1L)) {
      rows[[length(rows) + 1L]] <- data.table::data.table(
        covariate = paste0(nm, "__missing"),
        smd_before = weighted_smd_binary(miss, a, w0),
        smd_after = weighted_smd_binary(miss, a, wp)
      )
    }
  }
  for (nm in spec$categorical) {
    b <- categorical_level_smds(x[[nm]], a, w0, nm)
    p <- categorical_level_smds(x[[nm]], a, wp, nm)
    if (nrow(b) && nrow(p)) {
      data.table::setnames(b, "smd", "smd_before")
      data.table::setnames(p, "smd", "smd_after")
      rows[[length(rows) + 1L]] <- merge(b, p, by = "covariate", all = TRUE)
    }
  }
  out <- data.table::rbindlist(rows, fill = TRUE)
  out[, `:=`(abs_before = abs(smd_before), abs_after = abs(smd_after))]
  data.table::setorder(out, -abs_after)
  out
}

balance_summary <- function(tab) {
  data.table::data.table(
    statistic = c("max_abs_smd_before", "median_abs_smd_after", "max_abs_smd_after"),
    value = c(
      max(tab$abs_before, na.rm = TRUE),
      stats::median(tab$abs_after, na.rm = TRUE),
      max(tab$abs_after, na.rm = TRUE)
    )
  )
}
