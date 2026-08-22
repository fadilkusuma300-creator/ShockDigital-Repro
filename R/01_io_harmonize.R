# Input/output and harmonization helpers for World Bank Enterprise Survey files.
# The code resolves canonical variables through config/variable_map.yml and
# stops on ambiguous mappings rather than silently choosing a column.

source(file.path(project_root(), "R", "00_utils.R"))
source(file.path(project_root(), "R", "01a_credit_constraints.R"))

read_any_survey_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "dta") {
    return(data.table::as.data.table(haven::read_dta(path)))
  }
  if (ext == "sav") {
    return(data.table::as.data.table(haven::read_sav(path)))
  }
  if (ext %in% c("csv", "txt")) {
    return(data.table::fread(path, na.strings = c("", "NA", "N/A", ".")))
  }
  if (ext == "rds") {
    return(data.table::as.data.table(readRDS(path)))
  }
  stop("Unsupported survey file type: ", path)
}

list_survey_files <- function(directory) {
  if (!dir.exists(directory)) stop("Missing input directory: ", directory)
  files <- list.files(directory, pattern = "\\.(dta|sav|csv|txt|rds)$",
                      full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No supported survey files found in: ", directory)
  sort(files)
}

column_labels <- function(dt) {
  vapply(names(dt), function(nm) {
    lab <- attr(dt[[nm]], "label", exact = TRUE)
    if (is.null(lab)) "" else as.character(lab)
  }, character(1))
}

resolve_variable <- function(dt, spec, canonical_name, required = TRUE) {
  # 1) exact candidate-name matching, case-insensitive
  nms <- names(dt)
  lower <- tolower(nms)
  candidates <- tolower(unlist(spec$candidates %||% character()))
  hits <- which(lower %in% candidates)

  if (length(hits) == 1L) return(nms[hits])
  if (length(hits) > 1L) {
    # Prefer the first candidate by explicit configuration order.
    ordered_hits <- unlist(lapply(candidates, function(x) which(lower == x)))
    ordered_hits <- unique(ordered_hits)
    if (length(ordered_hits)) return(nms[ordered_hits[1]])
  }

  # 2) variable-label regular expressions
  regexes <- unlist(spec$label_regex %||% character())
  if (length(regexes)) {
    labs <- column_labels(dt)
    matches <- integer()
    for (rx in regexes) {
      matches <- union(matches, grep(rx, labs, ignore.case = TRUE, perl = TRUE))
    }
    if (length(matches) == 1L) return(nms[matches])
    if (length(matches) > 1L) {
      stop(
        "Ambiguous label mapping for '", canonical_name, "': ",
        paste(nms[matches], collapse = ", "),
        ". Set an exact candidate in config/variable_map.yml."
      )
    }
  }

  if (required) {
    stop(
      "Cannot resolve required variable '", canonical_name,
      "'. Add the source column name to config/variable_map.yml."
    )
  }
  NA_character_
}

resolve_map <- function(dt, section, required_names = names(section)) {
  out <- setNames(rep(NA_character_, length(section)), names(section))
  for (nm in names(section)) {
    out[[nm]] <- resolve_variable(
      dt, section[[nm]], nm,
      required = nm %in% required_names
    )
  }
  out
}

extract_canonical <- function(dt, mapping, missing_codes, source_file = NA_character_) {
  out <- data.table::data.table(row_id_source = seq_len(nrow(dt)))
  for (canonical in names(mapping)) {
    src <- mapping[[canonical]]
    if (is.na(src) || !nzchar(src)) {
      out[[canonical]] <- NA
    } else {
      out[[canonical]] <- missing_to_na(dt[[src]], missing_codes)
    }
  }
  out[, source_file := source_file]
  out
}

infer_round_from_filename <- function(path) {
  bn <- basename(path)
  patterns <- c(
    "(?:round|rnd|wave|r)[ _.-]*([0-9]+)",
    "(?:covid|followup)[ _.-]*([0-9]+)"
  )
  for (rx in patterns) {
    m <- regexec(rx, bn, ignore.case = TRUE, perl = TRUE)
    z <- regmatches(bn, m)[[1]]
    if (length(z) >= 2L) {
      candidate <- suppressWarnings(as.integer(z[2]))
      # Follow-up wave numbers are small positive integers. Capping the
      # filename-derived candidate prevents a calendar year such as 2020 from
      # being misread as a survey round.
      if (is.finite(candidate) && candidate >= 1L && candidate <= 20L) return(candidate)
    }
  }
  NA_integer_
}

harmonize_baseline <- function(root = project_root(), cfg = load_analysis_config(root),
                               vmap = load_variable_map(root)) {
  files <- list_survey_files(file.path(root, cfg$paths$baseline_dir))
  required <- c(
    "firm_id", "survey_weight", "baseline_employment", "establishment_year",
    "manager_experience", "website", "foreign_ownership_share",
    "direct_export_share", "indirect_export_share", "industry", "economy",
    "interview_year"
  )
  # Innovation is a direct baseline item. Credit constraints are handled
  # separately: a validated direct indicator is accepted, otherwise the project
  # derives the World Bank fin23/fin24/fin25 classification from K16/K17/K20a1
  # and K3/K5 external-finance components.
  if (isTRUE(cfg$runtime$strict_mapping)) {
    required <- c(required, "innovation")
  }

  pieces <- lapply(files, function(path) {
    raw <- read_any_survey_file(path)
    mapping <- resolve_map(raw, vmap$baseline, required_names = required)
    can <- extract_canonical(raw, mapping, vmap$coding$missing_codes, basename(path))
    attr(can, "mapping") <- mapping
    can
  })
  dt <- data.table::rbindlist(pieces, fill = TRUE, use.names = TRUE)

  # Enforce stable string identifiers before joining across survey vintages.
  dt[, firm_id := trimws(safe_character(firm_id))]
  dt[, economy := trimws(safe_character(economy))]
  dt[, industry := trimws(safe_character(industry))]

  numeric_cols <- c(
    "survey_weight", "baseline_employment", "establishment_year",
    "manager_experience", "foreign_ownership_share", "direct_export_share",
    "indirect_export_share", "interview_year", "credit_constraints", "innovation",
    "loan_applied", "loan_noapply_reason", "loan_application_outcome",
    "wc_bank_combined", "wc_bank_private", "wc_bank_state", "wc_nonbank",
    "wc_supplier_credit", "wc_other_external", "fa_bank_combined",
    "fa_bank_private", "fa_bank_state", "fa_nonbank", "fa_supplier_credit",
    "fa_other_external", "fa_equity"
  )
  for (nm in intersect(numeric_cols, names(dt))) dt[, (nm) := safe_numeric(get(nm))]

  # Construct the credit-constraint covariate from the official World Bank
  # classification whenever a validated direct indicator is not present.
  dt <- construct_credit_constraints(dt, strict = isTRUE(cfg$runtime$strict_mapping))

  # Common binary constructions used by the paper's baseline adjustment set.
  dt[, foreign_owned := data.table::fcase(
    is.finite(foreign_ownership_share), as.integer(foreign_ownership_share > 0),
    default = NA_integer_
  )]
  dt[, exporter := data.table::fcase(
    (is.finite(direct_export_share) & direct_export_share > 0) |
      (is.finite(indirect_export_share) & indirect_export_share > 0), 1L,
    is.finite(direct_export_share) & is.finite(indirect_export_share) &
      direct_export_share == 0 & indirect_export_share == 0, 0L,
    default = NA_integer_
  )]
  dt[, website := safe_numeric(website)]

  # Establishment age is anchored to the pre-pandemic survey year.
  dt[, firm_age := ifelse(
    is.finite(interview_year) & is.finite(establishment_year),
    pmax(0, interview_year - establishment_year), NA_real_
  )]

  # If duplicate baseline rows exist for an establishment, use the latest
  # pre-pandemic record because it is temporally closest to the follow-up. Firm
  # identifiers are not assumed to be globally unique across economies.
  data.table::setorder(dt, economy, firm_id, -interview_year)
  dt <- unique(dt, by = c("economy", "firm_id"))
  dt
}

harmonize_followup <- function(root = project_root(), cfg = load_analysis_config(root),
                               vmap = load_variable_map(root)) {
  files <- list_survey_files(file.path(root, cfg$paths$followup_dir))
  required <- c(
    "firm_id", "survey_weight", "economy", "interview_year", "interview_month",
    "operating_status", "sales_direction", "sales_increase_pct",
    "sales_decrease_pct", "demand_state", "input_supply_state",
    "employment_current"
  )

  pieces <- lapply(files, function(path) {
    raw <- read_any_survey_file(path)
    # Round, cash-flow variants, and online variants may legitimately be absent
    # in a particular questionnaire wave, so they are resolved as optional.
    mapping <- resolve_map(raw, vmap$followup, required_names = required)
    can <- extract_canonical(raw, mapping, vmap$coding$missing_codes, basename(path))
    can[, round_from_file := infer_round_from_filename(path)]
    attr(can, "mapping") <- mapping
    can
  })
  dt <- data.table::rbindlist(pieces, fill = TRUE, use.names = TRUE)

  dt[, firm_id := trimws(safe_character(firm_id))]
  dt[, economy := trimws(safe_character(economy))]
  numeric_cols <- setdiff(names(vmap$followup), c("firm_id", "economy", "strata"))
  for (nm in intersect(numeric_cols, names(dt))) dt[, (nm) := safe_numeric(get(nm))]

  # Prefer an explicit round variable. If it is unavailable, use a round number
  # embedded in the filename. A final date-based fallback is applied below.
  dt[, round_final := as.integer(round)]
  dt[!is.finite(round_final), round_final := as.integer(round_from_file)]

  # The public follow-up questionnaires use COVc4a in Round 1 and round-specific
  # variants such as COV2c4a/COV3c4a in later rounds. The same pattern holds
  # for the cash-flow item.
  dt[, online_adjust_raw := data.table::fifelse(
    round_final == 1L,
    safe_numeric(online_adjust_round1),
    safe_numeric(online_adjust_later),
    na = NA_real_
  )]
  dt[is.na(online_adjust_raw), online_adjust_raw := data.table::fcoalesce(
    safe_numeric(online_adjust_later), safe_numeric(online_adjust_round1)
  )]

  dt[, cash_flow_raw := data.table::fifelse(
    round_final == 1L,
    safe_numeric(cash_flow_state_round1),
    safe_numeric(cash_flow_state_later),
    na = NA_real_
  )]
  dt[is.na(cash_flow_raw), cash_flow_raw := data.table::fcoalesce(
    safe_numeric(cash_flow_state_later), safe_numeric(cash_flow_state_round1)
  )]

  dt[, interview_date_index := 12 * safe_numeric(interview_year) + safe_numeric(interview_month)]

  # Date-based round inference is only used if the supplied files contain no
  # usable round field or filename marker for a row. It ranks survey dates
  # within each economy and preserves ties.
  missing_round <- which(!is.finite(dt$round_final))
  if (length(missing_round)) {
    dt[, inferred_date_round := data.table::frank(
      interview_date_index, ties.method = "dense", na.last = "keep"
    ), by = economy]
    dt[missing_round, round_final := as.integer(inferred_date_round)]
  }

  code <- vmap$coding
  dt[, online_adjust := data.table::fcase(
    online_adjust_raw == code$yes, 1L,
    online_adjust_raw == code$no, 0L,
    default = NA_integer_
  )]

  # Signed sales change: increases positive, unchanged zero, decreases negative.
  dt[, sales_change := data.table::fcase(
    sales_direction == code$increase, abs(safe_numeric(sales_increase_pct)),
    sales_direction == code$same, 0,
    sales_direction == code$decrease, -abs(safe_numeric(sales_decrease_pct)),
    default = NA_real_
  )]

  dt[, operating_open := data.table::fcase(
    operating_status == code$open, 1L,
    operating_status %in% c(code$temporarily_closed, code$permanently_closed), 0L,
    default = NA_integer_
  )]

  dt[, employment_current := safe_numeric(employment_current)]
  dt[, demand_state := safe_numeric(demand_state)]
  dt[, input_supply_state := safe_numeric(input_supply_state)]
  dt[, cash_flow_raw := safe_numeric(cash_flow_raw)]

  # Remove exact duplicate firm-round records. Non-identical duplicates stop the
  # analysis because resolving them requires source-specific adjudication.
  key <- c("economy", "firm_id", "round_final")
  dup <- dt[, .N, by = key][N > 1L]
  if (nrow(dup)) {
    stop(
      "Duplicate firm-round records detected after harmonization. Resolve them ",
      "in the raw-data layer before analysis. Example: ",
      paste(unlist(dup[1, ..key]), collapse = " / ")
    )
  }
  data.table::setorder(dt, economy, firm_id, round_final, interview_date_index)
  dt
}
