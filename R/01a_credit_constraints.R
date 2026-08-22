# World Bank credit-constraint construction used for the baseline covariate.
#
# The classification follows the World Bank Enterprise Surveys finance
# indicators fin23/fin24/fin25 (Islam & Rodriguez Meza, 2023). The binary
# analysis covariate equals 1 for partially OR fully credit constrained and
# 0 for credit unconstrained. When the source survey does not provide enough
# finance items to classify a firm, the result is left missing rather than
# guessed.

source(file.path(project_root(), "R", "00_utils.R"))

row_any_positive <- function(dt, cols) {
  cols <- intersect(cols, names(dt))
  if (!length(cols)) return(rep(NA_integer_, nrow(dt)))
  m <- do.call(cbind, lapply(cols, function(nm) safe_numeric(dt[[nm]])))
  if (is.null(dim(m))) m <- matrix(m, ncol = 1L)
  any_obs <- rowSums(is.finite(m)) > 0L
  any_pos <- rowSums(is.finite(m) & m > 0, na.rm = TRUE) > 0L
  out <- rep(NA_integer_, nrow(dt))
  out[any_obs & !any_pos] <- 0L
  out[any_pos] <- 1L
  out
}

normalize_direct_credit_indicator <- function(x) {
  z <- safe_numeric(x)
  vals <- sort(unique(z[is.finite(z)]))
  if (!length(vals)) return(rep(NA_integer_, length(z)))
  if (all(vals %in% c(0, 1))) return(as.integer(z))
  if (all(vals %in% c(1, 2))) {
    return(data.table::fcase(z == 1, 1L, z == 2, 0L, default = NA_integer_))
  }
  stop(
    "A direct credit-constraint variable was found, but its coding is not 0/1 or 1/2. ",
    "Map the raw finance components and let the project derive the World Bank classification instead."
  )
}

derive_world_bank_credit_constraints <- function(dt) {
  # Current and legacy WBES files can split bank finance into private/state
  # columns or store it in a combined column. Every available external-finance
  # component is therefore inspected, and no missing component is silently
  # treated as zero.
  non_equity_cols <- c(
    "wc_bank_combined", "wc_bank_private", "wc_bank_state", "wc_nonbank",
    "wc_supplier_credit", "wc_other_external", "fa_bank_combined",
    "fa_bank_private", "fa_bank_state", "fa_nonbank", "fa_supplier_credit",
    "fa_other_external"
  )
  ext <- row_any_positive(dt, non_equity_cols)
  equity <- row_any_positive(dt, "fa_equity")

  applied <- safe_numeric(dt$loan_applied)
  reason <- safe_numeric(dt$loan_noapply_reason)
  outcome <- safe_numeric(dt$loan_application_outcome)

  # K16: 1 yes, 2 no. K17: 1 no need; 2-7 other reasons.
  # K20a1: 1 fully approved, 2 partially approved, 3 rejected.
  unconstrained <-
    (applied == 1 & outcome == 1) |
    (applied == 2 & reason == 1)

  partially <-
    (applied == 1 & outcome == 2) |
    (applied == 1 & outcome == 3 & ext == 1) |
    (applied == 2 & reason %in% 2:7 & ext == 1)

  fully <-
    (applied == 1 & outcome == 3 & ext == 0) |
    (applied == 2 & reason %in% 2:7 & ext == 0)

  cls <- rep(NA_character_, nrow(dt))
  cls[unconstrained %in% TRUE] <- "unconstrained"
  cls[partially %in% TRUE] <- "partially_constrained"
  cls[fully %in% TRUE] <- "fully_constrained"

  # Equity is retained as a classification field because the World Bank definition
  # explicitly distinguishes equity from the non-equity external-finance test.
  data.table::data.table(
    credit_constraint_class = cls,
    credit_constraints = data.table::fcase(
      cls %in% c("partially_constrained", "fully_constrained"), 1L,
      cls == "unconstrained", 0L,
      default = NA_integer_
    ),
    credit_external_finance_non_equity = ext,
    credit_equity_finance = equity
  )
}

construct_credit_constraints <- function(dt, strict = TRUE) {
  x <- data.table::copy(dt)

  direct <- rep(NA_integer_, nrow(x))
  if ("credit_constraints" %in% names(x) && any(is.finite(safe_numeric(x$credit_constraints)))) {
    direct <- normalize_direct_credit_indicator(x$credit_constraints)
  }

  required_components <- c("loan_applied", "loan_noapply_reason", "loan_application_outcome")
  raw_core_available <- all(required_components %in% names(x)) &&
    any(vapply(required_components, function(nm) any(is.finite(safe_numeric(x[[nm]]))), logical(1)))

  derived <- data.table::data.table(
    credit_constraint_class = rep(NA_character_, nrow(x)),
    credit_constraints = rep(NA_integer_, nrow(x)),
    credit_external_finance_non_equity = rep(NA_integer_, nrow(x)),
    credit_equity_finance = rep(NA_integer_, nrow(x))
  )
  if (raw_core_available) derived <- derive_world_bank_credit_constraints(x)

  # Use a validated direct indicator where present, then fill remaining rows
  # from the official raw-item construction. This supports projects that combine
  # survey vintages with different data layouts without sacrificing information.
  final <- direct
  fill <- is.na(final) & !is.na(derived$credit_constraints)
  final[fill] <- derived$credit_constraints[fill]

  cls <- rep(NA_character_, nrow(x))
  cls[!is.na(direct) & direct == 1L] <- "constrained_direct"
  cls[!is.na(direct) & direct == 0L] <- "unconstrained_direct"
  cls[is.na(direct) & !is.na(derived$credit_constraint_class)] <-
    derived$credit_constraint_class[is.na(direct) & !is.na(derived$credit_constraint_class)]

  source <- rep("unavailable", nrow(x))
  source[!is.na(direct)] <- "direct_source_indicator"
  source[is.na(direct) & !is.na(derived$credit_constraints)] <-
    "derived_world_bank_fin23_fin24_fin25"

  x[, credit_constraints := final]
  x[, credit_constraint_class := cls]
  x[, credit_constraint_source := source]
  x[, credit_external_finance_non_equity := derived$credit_external_finance_non_equity]
  x[, credit_equity_finance := derived$credit_equity_finance]

  if (isTRUE(strict) && !any(!is.na(x$credit_constraints))) {
    stop(
      "Cannot construct the credit-constraint covariate. Supply either a validated direct indicator or ",
      "the World Bank finance items K16, K17, K20a1 plus available K3/K5 external-finance components."
    )
  }
  x
}

