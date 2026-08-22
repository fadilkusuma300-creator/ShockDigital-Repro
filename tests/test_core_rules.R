# Lightweight smoke tests using only base R + project helpers.
# Run from the project root with: Rscript tests/test_core_rules.R
source("R/00_utils.R")
source("R/01_io_harmonize.R")
source("R/01a_credit_constraints.R")

stopifnot(infer_round_from_filename("Zambia-Round-3-COVID-follow-up.dta") == 3L)
stopifnot(infer_round_from_filename("survey_R2_file.dta") == 2L)
stopifnot(is.na(infer_round_from_filename("survey_2020_file.dta")))

stopifnot(identical(clip(c(-1, 0.5, 2), 0, 1), c(0, 0.5, 1)))
stopifnot(abs(weighted_mean(c(0, 2), c(1, 3)) - 1.5) < 1e-12)

cat("Core smoke tests passed.\n")


# World Bank credit-constraint classification smoke test.
syn <- data.table::data.table(
  loan_applied = c(1, 1, 1, 2, 2),
  loan_noapply_reason = c(NA, NA, NA, 1, 3),
  loan_application_outcome = c(1, 2, 3, NA, NA),
  wc_bank_combined = c(0, 0, 0, 0, 10),
  wc_nonbank = 0, wc_supplier_credit = 0, wc_other_external = 0,
  fa_bank_combined = 0, fa_nonbank = 0, fa_supplier_credit = 0,
  fa_other_external = 0, fa_equity = c(0, 0, 1, 0, 0)
)
cc <- derive_world_bank_credit_constraints(syn)
# A rejected applicant with equity finance but no non-equity external finance
# remains fully constrained under the World Bank definition.
stopifnot(identical(cc$credit_constraints, c(0L, 1L, 1L, 0L, 1L)))
stopifnot(identical(cc$credit_constraint_class,
  c("unconstrained", "partially_constrained", "fully_constrained", "unconstrained", "partially_constrained")))

cat("Credit-constraint smoke tests passed.\n")
