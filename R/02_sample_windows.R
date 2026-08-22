# Construction of the strict adjacent three-wave longitudinal sample.

source(file.path(project_root(), "R", "00_utils.R"))

prepare_survey_weights <- function(dt) {
  # Normalize weights within economy x round and winsorize at the 1st/99th percentiles.
  out <- data.table::copy(dt)
  out[, survey_weight_raw := safe_numeric(survey_weight)]
  out[, survey_weight_normalized := normalize_weights(survey_weight_raw),
      by = .(economy, round_final)]
  out[, analysis_weight := {
    w <- survey_weight_normalized
    if (all(!is.finite(w))) rep(1, .N) else winsorize_weight(w, c(0.01, 0.99))
  }, by = .(economy, round_final)]
  out
}

merge_baseline_followup <- function(baseline, followup, cfg, apply_sme_filter = TRUE) {
  b <- data.table::copy(baseline)
  f <- data.table::copy(followup)

  if (isTRUE(apply_sme_filter)) {
    b <- b[
      is.finite(baseline_employment) &
        baseline_employment >= cfg$sample$min_baseline_employment &
        baseline_employment <= cfg$sample$max_baseline_employment
    ]
  }

  # Firm identifiers can repeat across economies. Join on both whenever economy
  # is available in the baseline file.
  join_keys <- c("economy", "firm_id")
  keep_b <- unique(c(
    join_keys, "baseline_employment", "firm_age", "manager_experience",
    "website", "credit_constraints", "innovation", "foreign_owned", "exporter",
    "industry", "strata", "interview_year"
  ))
  keep_b <- intersect(keep_b, names(b))
  data.table::setnames(b, "interview_year", "baseline_interview_year", skip_absent = TRUE)
  keep_b[keep_b == "interview_year"] <- "baseline_interview_year"

  merged <- merge(f, b[, ..keep_b], by = join_keys, all = FALSE, suffixes = c("", "_baseline"))
  merged
}

add_relative_interview_month <- function(dt) {
  out <- data.table::copy(dt)
  # Relative time is economy-specific and measured from that economy's first
  # observed follow-up interview month.
  out[, first_economy_month := suppressWarnings(min(interview_date_index, na.rm = TRUE)),
      by = economy]
  out[!is.finite(first_economy_month), first_economy_month := NA_real_]
  out[, relative_interview_month := interview_date_index - first_economy_month]
  out[, first_economy_month := NULL]
  out
}

construct_adjacent_windows <- function(panel, cfg) {
  x <- data.table::copy(panel)
  data.table::setorder(x, economy, firm_id, round_final, interview_date_index)

  # The treatment is only eligible when no earlier observed follow-up round has
  # reported starting/increasing online activity. Prior missing treatment values
  # are not treated as "no"; an unknown history makes eligibility indeterminate.
  x[, prior_online_any := {
    v <- online_adjust
    out <- rep(NA_integer_, .N)
    if (.N >= 1L) out[1L] <- 0L
    if (.N > 1L) {
      for (j in 2:.N) {
        prev <- v[seq_len(j - 1L)]
        out[j] <- if (any(is.na(prev))) NA_integer_ else as.integer(any(prev == 1L))
      }
    }
    out
  }, by = .(economy, firm_id)]

  # Create all adjacent triples within each establishment. Adjacency is defined
  # by the ordered observed rounds; no interpolation across a missing round.
  x[, next_round := data.table::shift(round_final, type = "lead"), by = .(economy, firm_id)]
  x[, next2_round := data.table::shift(round_final, n = 2L, type = "lead"), by = .(economy, firm_id)]
  x[, next_online := data.table::shift(online_adjust, type = "lead"), by = .(economy, firm_id)]
  x[, next2_sales := data.table::shift(sales_change, n = 2L, type = "lead"), by = .(economy, firm_id)]
  x[, next2_open := data.table::shift(operating_open, n = 2L, type = "lead"), by = .(economy, firm_id)]
  x[, next2_employment := data.table::shift(employment_current, n = 2L, type = "lead"), by = .(economy, firm_id)]
  x[, next2_input := data.table::shift(input_supply_state, n = 2L, type = "lead"), by = .(economy, firm_id)]
  x[, t1_rel_month := data.table::shift(relative_interview_month, type = "lead"), by = .(economy, firm_id)]
  x[, t2_rel_month := data.table::shift(relative_interview_month, n = 2L, type = "lead"), by = .(economy, firm_id)]
  x[, t1_round := data.table::shift(round_final, type = "lead"), by = .(economy, firm_id)]
  x[, t2_round := data.table::shift(round_final, n = 2L, type = "lead"), by = .(economy, firm_id)]

  # Preserve t0 state variables for adjustment and outcome construction.
  x[, shock_sales_decline := pmax(-sales_change, 0)]
  x[, t0_sales_change := sales_change]
  x[, t0_relative_interview_month := relative_interview_month]
  x[, t0_employment := employment_current]
  x[, t0_input_state := input_supply_state]
  x[, t0_demand_state := demand_state]
  x[, t0_cash_flow_state := cash_flow_raw]
  x[, t0_operating_status := operating_status]
  x[, t0_open := operating_open]

  # Strictly adjacent survey round numbers. This excludes a 1->3->4 sequence.
  x[, consecutive_triple :=
      is.finite(round_final) & is.finite(next_round) & is.finite(next2_round) &
      next_round == round_final + 1L & next2_round == round_final + 2L]

  # At t0 the establishment must not already have reported online expansion.
  # At t1, treatment is the first observed start/increase; control is no change.
  x[, history_clear_at_t0 :=
      !is.na(prior_online_any) & prior_online_any == 0L & online_adjust == 0L]

  cand <- x[
    consecutive_triple == TRUE &
      history_clear_at_t0 == TRUE &
      !is.na(next_online) &
      is.finite(shock_sales_decline)
  ]
  cand[, treatment := as.integer(next_online == 1L)]

  # If multiple eligible windows remain, retain the earliest.
  data.table::setorder(cand, economy, firm_id, round_final)
  cand <- cand[, .SD[1L], by = .(economy, firm_id)]

  # Primary and secondary outcomes at t2.
  cand[, y_sales := next2_sales]
  cand[, y_continued_operation := next2_open]
  cand[, y_employment_retention := data.table::fcase(
    is.finite(next2_employment) & is.finite(t0_employment),
      as.integer(next2_employment >= t0_employment),
    default = NA_integer_
  )]

  # Input states are coded increase / same / decrease in the public instrument.
  # Primary recovery means the t2 input condition is no longer worse than t0;
  # the strict definition requires an increase at t2.
  code <- load_variable_map()$coding
  cand[, y_input_recovery := data.table::fcase(
    is.finite(t0_input_state) & is.finite(next2_input),
      as.integer(next2_input != code$decrease),
    default = NA_integer_
  )]
  cand[, y_input_increase := data.table::fcase(
    is.finite(next2_input), as.integer(next2_input == code$increase),
    default = NA_integer_
  )]

  # Alternative sales outcome: leave the severe-sales-decline state by t2.
  cand[, severe_shock := as.integer(shock_sales_decline >= cfg$sample$severe_shock_cutoff)]
  cand[, y_exit_severe_sales_decline := data.table::fcase(
    is.finite(next2_sales),
      as.integer(pmax(-next2_sales, 0) < cfg$sample$severe_shock_cutoff),
    default = NA_integer_
  )]

  # Timing covariates are anchored at t0 so they are strictly pre-treatment.
  cand[, survey_round := round_final]
  cand[, relative_interview_month := t0_relative_interview_month]
  cand[]
}

sample_flow_table <- function(baseline, followup, matched_all, merged_sme, windows, cfg) {
  data.table::data.table(
    stage = c(
      "followup_firm_round_records",
      "followup_unique_firms",
      "matched_baseline_unique_firms",
      "sme_longitudinal_pool",
      "strict_adjacent_three_wave_sample",
      "economies_in_strict_sample",
      "treated_first_online_expansion",
      "control_no_online_expansion",
      "severe_t0_sales_shock"
    ),
    observed = c(
      nrow(followup),
      data.table::uniqueN(followup, by = c("economy", "firm_id")),
      data.table::uniqueN(matched_all, by = c("economy", "firm_id")),
      data.table::uniqueN(merged_sme, by = c("economy", "firm_id")),
      nrow(windows),
      data.table::uniqueN(windows$economy),
      sum(windows$treatment == 1L, na.rm = TRUE),
      sum(windows$treatment == 0L, na.rm = TRUE),
      sum(windows$severe_shock == 1L, na.rm = TRUE)
    )
  )
}

# Detailed sample-flow summary for longitudinal eligibility.
detailed_sample_flow_table <- function(merged_sme, windows) {
  firm <- merged_sme[, {
    rr <- sort(unique(round_final[is.finite(round_final)]))
    has_consecutive <- FALSE
    if (length(rr) >= 3L) {
      has_consecutive <- any(vapply(rr, function(r) all(c(r, r + 1L, r + 2L) %in% rr), logical(1)))
    }
    list(n_distinct_rounds = length(rr), has_consecutive_three = has_consecutive)
  }, by = .(economy, firm_id)]

  data.table::data.table(
    stage = c(
      "sme_matched_to_any_followup",
      "sme_with_at_least_three_distinct_followup_rounds",
      "sme_with_at_least_one_consecutive_three_round_sequence",
      "strict_eligible_earliest_window"
    ),
    observed = c(
      nrow(firm),
      sum(firm$n_distinct_rounds >= 3L),
      sum(firm$has_consecutive_three),
      nrow(windows)
    )
  )
}
