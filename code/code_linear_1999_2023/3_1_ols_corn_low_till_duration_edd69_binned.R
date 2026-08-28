## =============================== PORTABLE PATHS ==============================
## Set ROOT_OVERRIDE to the folder that contains "data/" and leave everything else
## alone. If you leave it as "", the script locates the project automatically.
##
##
## Required layout under that folder:
##   data/Field_yield_tillage/processed_smooth_1999_2023/<crop>_yield_tillage_XX_1999_2023.csv
##   data/PRISM/PRISM_climate_1999_2023.csv
##   data/irrigation/field_irrigation_ever_only_irrigated.csv
##
## Automatic detection tries, in order: ROOT_OVERRIDE, the PROJECT_ROOT environment
## variable, the folder of this script (works when source()d in RStudio and when run
## with Rscript), then the working directory -- ascending up to 12 levels from each.
## ============================================================================
ROOT_OVERRIDE <- ""

library(data.table)
library(fixest)
library(ggplot2)

find_root_dir <- function() {
  ascend_paths <- function(path, max_up = 12L) {
    out <- character()
    cur <- normalizePath(path, winslash = "/", mustWork = FALSE)
    for (i in seq_len(max_up)) {
      if (nzchar(cur)) out <- c(out, cur)
      parent <- dirname(cur)
      if (identical(parent, cur)) break
      cur <- parent
    }
    unique(out)
  }

  ## the project is identified by its DATA, not by a code/ subfolder, so the
  ## scripts can sit anywhere relative to the data
  looks_like_root <- function(p) {
    dir.exists(file.path(p, "data", "Field_yield_tillage",
                         "processed_smooth_1999_2023")) &&
      file.exists(file.path(p, "data", "PRISM", "PRISM_climate_1999_2023.csv"))
  }

  candidates <- character()

  ## 1. explicit override
  if (exists("ROOT_OVERRIDE", inherits = TRUE)) {
    ov <- get("ROOT_OVERRIDE", inherits = TRUE)
    if (is.character(ov) && length(ov) == 1L && nzchar(ov)) {
      if (!looks_like_root(ov))
        stop("ROOT_OVERRIDE is set to '", ov, "' but that folder does not contain\n",
             "  data/Field_yield_tillage/processed_smooth_1999_2023/  and\n",
             "  data/PRISM/PRISM_climate_1999_2023.csv\n",
             "Fix the path or set it back to \"\" for automatic detection.")
      return(normalizePath(ov))
    }
  }

  ## 2. environment variable
  env <- Sys.getenv("PROJECT_ROOT", unset = "")
  if (nzchar(env)) candidates <- c(candidates, ascend_paths(env))

  ## 3. this script's own folder -- source() in RStudio populates sys.frames()$ofile
  frame_files <- character()
  for (fr in sys.frames()) {
    of <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(of)) frame_files <- c(frame_files, as.character(of))
  }
  ## RStudio "Run selection" leaves no ofile; fall back to the editor context
  if (!length(frame_files) && requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (is.character(p) && length(p) == 1L && nzchar(p)) frame_files <- p
  }
  ## Rscript / R CMD BATCH
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) frame_files <- c(frame_files, sub("^--file=", "", fa[1]))
  for (f in unique(frame_files))
    candidates <- c(candidates, ascend_paths(dirname(f)))

  ## 4. working directory
  candidates <- unique(c(candidates, ascend_paths(getwd())))

  for (cand in candidates) if (looks_like_root(cand)) return(normalizePath(cand))

  stop("Could not locate the project root.\n",
       "Set ROOT_OVERRIDE at the top of this script to the folder containing data/.\n",
       "Searched ", length(candidates), " candidate folder(s), starting from:\n  ",
       paste(utils::head(candidates, 4), collapse = "\n  "))
}

sample_counties_by_state <- function(x, n_rows, seed) {
  if (is.na(n_rows) || n_rows <= 0 || n_rows >= nrow(x)) return(copy(x))

  state_targets <- x[, .(rows = .N), by = .(state = FIPS %/% 1000L)]
  state_targets[, target_rows := pmax(1L, as.integer(round(rows / sum(rows) * n_rows)))]
  sampled_parts <- vector("list", nrow(state_targets))

  for (i in seq_len(nrow(state_targets))) {
    state_id <- state_targets$state[i]
    state_target <- state_targets$target_rows[i]
    county_rows <- x[FIPS %/% 1000L == state_id, .(rows = .N), by = FIPS]
    set.seed(seed + i)
    county_rows <- county_rows[sample.int(.N)]
    county_rows[, cumulative_rows := cumsum(rows)]
    keep_fips <- county_rows[cumulative_rows <= state_target, FIPS]
    if (length(keep_fips) == 0) keep_fips <- county_rows[1, FIPS]
    sampled_parts[[i]] <- x[FIPS %in% keep_fips]
  }

  rbindlist(sampled_parts, use.names = TRUE)
}

make_tercile_bins <- function(x) {
  ranks <- frank(x, ties.method = "average", na.last = "keep")
  rank_pct <- ranks / sum(!is.na(x))
  bins <- cut(
    rank_pct,
    breaks = c(0, 1 / 4, 3 / 4, Inf),
    include.lowest = TRUE,
    labels = c("low", "middle", "high")
  )
  factor(as.character(bins), levels = c("low", "middle", "high"))
}

extract_term_table <- function(model, model_name, weather_var) {
  tab <- as.data.table(coeftable(model), keep.rownames = "term")
  setnames(
    tab,
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "std_error", "t_value", "p_value")
  )
  tab[, `:=`(model = model_name, weather_var = weather_var)]
  tab[]
}

theme_binning <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
      axis.title = element_text(face = "bold", color = "grey20"),
      axis.text = element_text(color = "grey25"),
      strip.text = element_text(face = "bold", color = "grey15"),
      strip.background = element_rect(fill = "grey94", color = NA),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", color = "grey10", size = base_size + 3),
      plot.subtitle = element_text(color = "grey35", size = base_size),
      plot.caption = element_text(color = "grey45", size = base_size - 2),
      plot.margin = margin(10, 14, 10, 10)
    )
}

estimate_binning_model <- function(dt, weather_var, weather_label, cluster_var, additive_weather_controls) {
  work <- copy(dt[is.finite(get(weather_var))])
  sample_bounds <- work[, .(
    weather_p05 = as.numeric(quantile(get(weather_var), 0.05, na.rm = TRUE)),
    weather_p95 = as.numeric(quantile(get(weather_var), 0.95, na.rm = TRUE))
  )]

  work[, weather_bin := make_tercile_bins(get(weather_var))]
  work <- work[!is.na(weather_bin)]

  bin_stats <- work[, .(
    observations = .N,
    fields = uniqueN(OBJECTID),
    counties = uniqueN(FIPS),
    control_observations = sum(low_till_duration == 0),
    low_till_observations = sum(low_till_duration > 0),
    mean_low_till_duration = as.numeric(mean(low_till_duration[low_till_duration > 0], na.rm = TRUE)),
    median_low_till_duration = as.numeric(median(low_till_duration[low_till_duration > 0], na.rm = TRUE)),
    weather_min = min(get(weather_var), na.rm = TRUE),
    weather_p05 = as.numeric(quantile(get(weather_var), 0.05, na.rm = TRUE)),
    weather_q25 = as.numeric(quantile(get(weather_var), 0.25, na.rm = TRUE)),
    weather_median = median(get(weather_var), na.rm = TRUE),
    weather_mean = mean(get(weather_var), na.rm = TRUE),
    weather_q75 = as.numeric(quantile(get(weather_var), 0.75, na.rm = TRUE)),
    weather_p95 = as.numeric(quantile(get(weather_var), 0.95, na.rm = TRUE)),
    weather_max = max(get(weather_var), na.rm = TRUE)
  ), by = weather_bin]

  work <- merge(
    work,
    bin_stats[, .(weather_bin, weather_bin_median = weather_median)],
    by = "weather_bin",
    all.x = TRUE,
    sort = FALSE
  )
  work[, weather_centered := get(weather_var) - weather_bin_median]

  bin_levels <- c("low", "middle", "high")
  active_bins <- setdiff(bin_levels, "low")

  for (b in bin_levels) {
    in_bin <- work$weather_bin == b
    if (b %in% active_bins) {
      work[, (paste0("bin_", b)) := as.integer(in_bin)]
    }
    work[, (paste0("weather_c_", b)) := as.integer(in_bin) * weather_centered]
    work[, (paste0("duration_years_", b)) := as.integer(in_bin) * low_till_duration]
    work[, (paste0("duration_years_x_weather_c_", b)) := low_till_duration * as.integer(in_bin) * weather_centered]
  }

  bin_intercept_terms <- paste0("bin_", active_bins)
  duration_terms <- paste0("duration_years_", bin_levels)
  weather_terms <- paste0("weather_c_", bin_levels)
  duration_slope_terms <- paste0("duration_years_x_weather_c_", bin_levels)
  rhs_terms <- c(bin_intercept_terms, duration_terms, weather_terms, duration_slope_terms, additive_weather_controls)

  model <- feols(
    as.formula(paste0("log_y ~ ", paste(rhs_terms, collapse = " + "), " | OBJECTID + year")),
    cluster = as.formula(paste0("~", cluster_var)),
    weights = ~size,
    lean = FALSE,
    mem.clean = TRUE,
    data = work
  )

  vc <- vcov(model)
  coefs <- coef(model)

  coef_and_se <- function(term) {
    if (!(term %in% names(coefs)) || !(term %in% rownames(vc)) || !(term %in% colnames(vc))) {
      return(c(estimate = NA_real_, std_error = NA_real_))
    }
    c(estimate = unname(coefs[term]), std_error = sqrt(vc[term, term]))
  }

  effects <- rbindlist(lapply(seq_along(bin_levels), function(i) {
    b <- bin_levels[i]
    term <- paste0("duration_years_", b)
    est <- coef_and_se(term)
    estimate <- est[["estimate"]]
    std_error <- est[["std_error"]]
    data.table(
      weather_var = weather_label,
      weather_source_var = weather_var,
      treatment = "low_till_duration_years",
      weather_bin = b,
      bin_order = i,
      estimate_log_points_per_year = estimate,
      std_error = std_error,
      ci_low_log_points_per_year = estimate - 1.96 * std_error,
      ci_high_log_points_per_year = estimate + 1.96 * std_error,
      pct_effect_per_year = (exp(estimate) - 1) * 100,
      pct_effect_low_per_year = (exp(estimate - 1.96 * std_error) - 1) * 100,
      pct_effect_high_per_year = (exp(estimate + 1.96 * std_error) - 1) * 100
    )
  }), use.names = TRUE)

  duration_weather_slopes <- rbindlist(lapply(seq_along(bin_levels), function(i) {
    b <- bin_levels[i]
    term <- paste0("duration_years_x_weather_c_", b)
    est <- coef_and_se(term)
    estimate <- est[["estimate"]]
    std_error <- est[["std_error"]]
    data.table(
      weather_var = weather_label,
      weather_source_var = weather_var,
      treatment = "low_till_duration_years",
      weather_bin = b,
      bin_order = i,
      slope_log_points_per_year_per_weather_unit = estimate,
      slope_std_error = std_error,
      slope_ci_low = estimate - 1.96 * std_error,
      slope_ci_high = estimate + 1.96 * std_error
    )
  }), use.names = TRUE)

  control_weather_slopes <- rbindlist(lapply(seq_along(bin_levels), function(i) {
    b <- bin_levels[i]
    term <- paste0("weather_c_", b)
    est <- coef_and_se(term)
    data.table(
      weather_var = weather_label,
      weather_source_var = weather_var,
      weather_bin = b,
      bin_order = i,
      control_weather_slope = est[["estimate"]],
      control_weather_slope_std_error = est[["std_error"]]
    )
  }), use.names = TRUE)

  additive_control_terms <- rbindlist(lapply(additive_weather_controls, function(term) {
    est <- coef_and_se(term)
    data.table(
      weather_var = weather_label,
      focal_weather_source_var = weather_var,
      additive_control = term,
      estimate = est[["estimate"]],
      std_error = est[["std_error"]],
      ci_low = est[["estimate"]] - 1.96 * est[["std_error"]],
      ci_high = est[["estimate"]] + 1.96 * est[["std_error"]]
    )
  }), use.names = TRUE, fill = TRUE)

  bin_stats[, weather_bin := as.character(weather_bin)]
  effects <- merge(effects, bin_stats, by = "weather_bin", all.x = TRUE, sort = FALSE)
  duration_weather_slopes <- merge(duration_weather_slopes, bin_stats, by = "weather_bin", all.x = TRUE, sort = FALSE)
  control_weather_slopes <- merge(control_weather_slopes, bin_stats, by = "weather_bin", all.x = TRUE, sort = FALSE)

  setorder(effects, bin_order)
  setorder(duration_weather_slopes, bin_order)
  setorder(control_weather_slopes, bin_order)

  list(
    model = model,
    effects = effects,
    duration_weather_slopes = duration_weather_slopes,
    control_weather_slopes = control_weather_slopes,
    additive_control_terms = additive_control_terms,
    terms = extract_term_table(model, paste0("linear_duration_binning_other_controls_", weather_label), weather_label),
    bin_stats = bin_stats[, `:=`(
      weather_var = weather_label,
      weather_source_var = weather_var,
      trimmed_sample_p05 = sample_bounds$weather_p05,
      trimmed_sample_p95 = sample_bounds$weather_p95
    )]
  )
}

make_hist_data <- function(dist, weather_id, y_floor, y_span, dist_hist_bins, weather_order, weather_labels) {
  empty <- data.table(
    weather_var = factor(character(), levels = weather_order),
    weather_label = factor(character(), levels = weather_labels[weather_order]),
    xmin = numeric(), xmax = numeric(), ymin = numeric(), ymax = numeric(), count = integer()
  )
  weather_id <- as.character(weather_id)
  if (is.null(dist) || nrow(dist) == 0) return(empty)
  d <- dist[weather_var == weather_id]
  if (nrow(d) == 0) return(empty)
  h_raw <- hist(d$weather_value, breaks = dist_hist_bins, plot = FALSE)
  h <- data.table(xmin = head(h_raw$breaks, -1), xmax = tail(h_raw$breaks, -1), count = h_raw$counts)
  h <- h[count > 0]
  if (nrow(h) == 0) return(empty)
  band_height <- 0.14 * y_span
  h[, `:=`(
    weather_var = weather_id,
    weather_label = factor(weather_labels[[weather_id]], levels = weather_labels[weather_order]),
    ymin = y_floor,
    ymax = y_floor + (count / max(count)) * band_height
  )]
  h[]
}

extract_fixed_effect_tables <- function(model_list) {
  field_fe_parts <- list()
  year_fe_parts <- list()
  summary_parts <- list()

  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    fe_list <- fixef(model)
    fe_names <- names(fe_list)

    summary_parts[[model_name]] <- data.table(
      model_name = model_name,
      fixed_effect = fe_names,
      n_levels = vapply(fe_list, length, integer(1))
    )

    if ("OBJECTID" %in% fe_names) {
      field_fe_parts[[model_name]] <- data.table(
        model_name = model_name,
        OBJECTID = names(fe_list[["OBJECTID"]]),
        field_fixed_effect = as.numeric(fe_list[["OBJECTID"]])
      )
    }

    if ("year" %in% fe_names) {
      year_fe_parts[[model_name]] <- data.table(
        model_name = model_name,
        year = as.integer(names(fe_list[["year"]])),
        year_fixed_effect = as.numeric(fe_list[["year"]])
      )
    }
  }

  list(
    field_fixed_effects = rbindlist(field_fe_parts, use.names = TRUE, fill = TRUE),
    year_fixed_effects = rbindlist(year_fe_parts, use.names = TRUE, fill = TRUE),
    fixed_effects_summary = rbindlist(summary_parts, use.names = TRUE, fill = TRUE)
  )
}

extract_vcov_tables <- function(model_list) {
  rbindlist(lapply(names(model_list), function(model_name) {
    vc <- vcov(model_list[[model_name]])
    out <- as.data.table(as.table(vc))
    setnames(out, c("term_1", "term_2", "covariance"))
    out[, model_name := model_name]
    out[]
  }), use.names = TRUE, fill = TRUE)
}

root_dir <- find_root_dir()
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "output_linear_1999_2023", "3_1_corn_low_till_duration_weather_binning_edd69")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(output_dir, "figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

crop <- "corn"
target_rows <- as.integer(Sys.getenv("TARGET_ROWS", unset = "50000000"))
sample_seed <- as.integer(Sys.getenv("SAMPLE_SEED", unset = "42"))
cluster_var <- Sys.getenv("CLUSTER_VAR", unset = "FIPS")
effect_scale <- Sys.getenv("EFFECT_SCALE", unset = "percent")
show_distribution <- as.logical(Sys.getenv("SHOW_DISTRIBUTION", unset = "TRUE"))
dist_hist_bins <- as.integer(Sys.getenv("DIST_HIST_BINS", unset = "35"))
exclude_years_raw <- Sys.getenv("EXCLUDE_YEARS", unset = "")
exclude_years <- if (nzchar(exclude_years_raw)) {
  as.integer(trimws(strsplit(exclude_years_raw, ",", fixed = TRUE)[[1]]))
} else {
  integer()
}
exclude_years <- sort(unique(exclude_years[is.finite(exclude_years)]))

message("Root directory: ", root_dir)
message("Crop: ", crop)
message("Target rows: ", target_rows)
message("Cluster variable: ", cluster_var)
message("Exclude years: ", if (length(exclude_years) == 0) "none" else paste(exclude_years, collapse = ", "))

field_dir <- file.path(data_dir, "Field_yield_tillage", "processed_smooth_1999_2023")
pattern <- paste0("^", crop, "_yield_tillage_[A-Z]{2}_1999_2023\\.csv$")
field_files <- list.files(field_dir, pattern = pattern, full.names = TRUE)
if (length(field_files) == 0) stop("No field files found for crop = ", crop)

field_cols <- c(
  "OBJECTID", "yield", "size", "year", "till", "GEOID",
  "till_1_count", "till_1_streak", "till_0_count", "till_0_streak"
)

dt <- rbindlist(lapply(field_files, fread, select = field_cols), use.names = TRUE)
setnames(dt, "GEOID", "FIPS")

dt[, `:=`(
  OBJECTID = as.character(OBJECTID),
  FIPS = as.integer(FIPS),
  year = as.integer(year),
  till = as.integer(till),
  till_1_count = as.integer(till_1_count),
  till_1_streak = as.integer(till_1_streak),
  till_0_count = as.integer(till_0_count),
  till_0_streak = as.integer(till_0_streak),
  size = as.numeric(size),
  yield_bu = as.numeric(yield) / 1000
)]

dt <- dt[is.finite(yield_bu) & yield_bu > 0 & is.finite(size) & size > 0]

omit.fields = c(unique(dt[year == 1999 &
                            till_1_count == 1]$OBJECTID), unique(dt[year == 2000 &
                                                                      till_1_count == 2]$OBJECTID))

dt <- dt[!OBJECTID %in% omit.fields]

dt[, log_y := log(yield_bu)]
if (length(exclude_years) > 0) dt <- dt[!(year %in% exclude_years)]

dt[, clean_control := till == 0L & till_1_count == 0L]
dt[, continuous_low_till := till == 1L & till_1_count == till_1_streak]
dt_duration <- dt[clean_control | continuous_low_till]
dt_duration[, low_till_duration := fifelse(clean_control, 0L, till_1_streak)]
dt_duration <- dt_duration[is.finite(low_till_duration) & low_till_duration >= 0]

climate_source_cols <- c("FIPS", "year", "GDD_4_5", "ppt_4_5", "GDD_6_9", "ppt_6_9", "EDD_6_9")
dt_climate <- fread(file.path(data_dir, "PRISM", "PRISM_climate_1999_2023.csv"), select = climate_source_cols)
dt_climate[, `:=`(
  FIPS = as.integer(FIPS),
  year = as.integer(year),
  GDD_4_5 = as.numeric(GDD_4_5),
  ppt_4_5 = as.numeric(ppt_4_5),
  GDD_6_9 = as.numeric(GDD_6_9),
  ppt_6_9 = as.numeric(ppt_6_9),
  EDD_6_9 = as.numeric(EDD_6_9)
)]

dt_duration <- merge(dt_duration, dt_climate, by = c("FIPS", "year"), all.x = TRUE)
weather_map <- c(EDD_6_9 = "EDD_6_9")
all_weather_controls <- c("GDD_4_5", "ppt_4_5", "GDD_6_9", "ppt_6_9", "EDD_6_9")

dt_duration <- dt_duration[
  is.finite(GDD_4_5) &
    is.finite(ppt_4_5) &
    is.finite(GDD_6_9) &
    is.finite(ppt_6_9) &
    is.finite(EDD_6_9)
]

dt_sample <- sample_counties_by_state(dt_duration, target_rows, sample_seed)
dt_sample[, state := as.factor(FIPS %/% 1000L)]

dt_irrigation <- fread(file.path(data_dir, "irrigation", "field_irrigation_ever_only_irrigated.csv"))
dt_sample <- dt_sample[!dt_sample$OBJECTID %in% dt_irrigation$OBJECTID, ]

if (!(cluster_var %in% names(dt_sample))) {
  stop("CLUSTER_VAR must be a column in the estimation data. Received: ", cluster_var)
}

sample_counts <- rbindlist(list(
  dt[, .(
    sample = "full_cleaned_panel",
    observations = .N,
    fields = uniqueN(OBJECTID),
    counties = uniqueN(FIPS),
    years = uniqueN(year),
    low_till_observations = NA_integer_,
    mean_low_till_duration = NA_real_
  )],
  dt_duration[, .(
    sample = "clean_controls_and_continuous_low_till_with_weather",
    observations = .N,
    fields = uniqueN(OBJECTID),
    counties = uniqueN(FIPS),
    years = uniqueN(year),
    low_till_observations = sum(low_till_duration > 0, na.rm = TRUE),
    mean_low_till_duration = mean(low_till_duration[low_till_duration > 0], na.rm = TRUE)
  )],
  dt_sample[, .(
    sample = "estimation_sample",
    observations = .N,
    fields = uniqueN(OBJECTID),
    counties = uniqueN(FIPS),
    years = uniqueN(year),
    low_till_observations = sum(low_till_duration > 0, na.rm = TRUE),
    mean_low_till_duration = mean(low_till_duration[low_till_duration > 0], na.rm = TRUE)
  )]
), use.names = TRUE)

results <- lapply(names(weather_map), function(label) {
  source_var <- unname(weather_map[label])
  additive_weather_controls <- setdiff(all_weather_controls, source_var)
  message(
    "Estimating focal weather-binning model for ", label, " (", source_var,
    ") with additive controls: ", paste(additive_weather_controls, collapse = ", ")
  )
  estimate_binning_model(dt_sample, source_var, label, cluster_var, additive_weather_controls)
})
names(results) <- names(weather_map)

effects <- rbindlist(lapply(results, `[[`, "effects"), use.names = TRUE, fill = TRUE)
duration_weather_slopes <- rbindlist(lapply(results, `[[`, "duration_weather_slopes"), use.names = TRUE, fill = TRUE)
control_weather_slopes <- rbindlist(lapply(results, `[[`, "control_weather_slopes"), use.names = TRUE, fill = TRUE)
additive_control_terms <- rbindlist(lapply(results, `[[`, "additive_control_terms"), use.names = TRUE, fill = TRUE)
term_tables <- rbindlist(lapply(results, `[[`, "terms"), use.names = TRUE, fill = TRUE)
bin_stats <- rbindlist(lapply(results, `[[`, "bin_stats"), use.names = TRUE, fill = TRUE)
model_list <- lapply(results, `[[`, "model")
fixed_effect_tables <- extract_fixed_effect_tables(model_list)
vcov_tables <- extract_vcov_tables(model_list)

fwrite(sample_counts, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_sample_counts.csv"))
fwrite(effects, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_effects.csv"))
fwrite(duration_weather_slopes, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_duration_weather_slopes.csv"))
fwrite(control_weather_slopes, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_control_weather_slopes.csv"))
fwrite(additive_control_terms, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_additive_control_terms.csv"))
fwrite(term_tables, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_terms.csv"))
fwrite(bin_stats, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_bin_stats.csv"))
fwrite(fixed_effect_tables$field_fixed_effects, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_field_fixed_effects.csv"))
fwrite(fixed_effect_tables$year_fixed_effects, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_year_fixed_effects.csv"))
fwrite(fixed_effect_tables$fixed_effects_summary, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_fixed_effects_summary.csv"))
fwrite(vcov_tables, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_vcov.csv"))
saveRDS(model_list, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_models.rds"))

weather_order <- c("EDD_6_9")
weather_labels <- c(EDD_6_9 = "EDD, Jun-Sep")
weather_units <- c(EDD_6_9 = "EDD")
exclude_years_note <- if (length(exclude_years) == 0) "Excluded years: none." else paste0("Excluded years: ", paste(exclude_years, collapse = ", "), ".")

plot_effects <- copy(effects)
plot_effects[, weather_var := factor(weather_var, levels = weather_order)]
plot_effects[, weather_label := factor(weather_labels[as.character(weather_var)], levels = weather_labels[weather_order])]
plot_effects[, weather_bin := factor(weather_bin, levels = c("low", "middle", "high"))]

if (effect_scale == "log") {
  plot_effects[, `:=`(
    estimate = estimate_log_points_per_year,
    ci_low = ci_low_log_points_per_year,
    ci_high = ci_high_log_points_per_year
  )]
  y_label <- "Effect per additional low-till year on corn yield"
  scale_suffix <- "log_points"
} else {
  plot_effects[, `:=`(
    estimate = pct_effect_per_year,
    ci_low = pct_effect_low_per_year,
    ci_high = pct_effect_high_per_year
  )]
  y_label <- "Effect per additional low-till year on corn yield (%)"
  scale_suffix <- "percent"
}

bin_ranges <- unique(plot_effects[, .(
  weather_var, weather_label, weather_bin, bin_order,
  weather_min, weather_median, weather_max,
  observations, control_observations, low_till_observations
)])
bin_ranges[, bin_label := paste0(weather_bin, "\nmedian=", round(weather_median, 1))]

distribution_data <- if (show_distribution) {
  vals <- dt_sample[is.finite(EDD_6_9), EDD_6_9]
  vals <- vals[vals >= quantile(vals, 0.025, na.rm = TRUE) & vals <= quantile(vals, 0.975, na.rm = TRUE)]
  data.table(weather_var = "EDD_6_9", weather_value = vals)
} else NULL

if (!is.null(distribution_data)) {
  distribution_data[, weather_var := factor(weather_var, levels = weather_order)]
}

make_single_weather_plot <- function(weather_id) {
  d <- plot_effects[weather_var == weather_id]
  r <- bin_ranges[weather_var == weather_id]
  y_floor <- min(d$ci_low, na.rm = TRUE)
  y_ceiling <- max(d$ci_high, na.rm = TRUE)
  y_span <- y_ceiling - y_floor
  hist_floor <- y_floor - 0.20 * y_span
  hist <- make_hist_data(distribution_data, weather_id, hist_floor, y_span, dist_hist_bins, weather_order, weather_labels)

  ggplot() +
    geom_rect(data = r, aes(xmin = weather_min, xmax = weather_max, ymin = -Inf, ymax = Inf, fill = weather_bin), alpha = 0.12, inherit.aes = FALSE) +
    geom_rect(data = hist, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey45", alpha = 0.30, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
    geom_errorbar(data = d, aes(x = weather_median, ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.85, color = "#276B3A") +
    geom_line(data = d, aes(x = weather_median, y = estimate, group = 1), linewidth = 0.7, color = "#276B3A") +
    geom_point(data = d, aes(x = weather_median, y = estimate), size = 2.8, stroke = 0.7, color = "#276B3A") +
    geom_text(data = r, aes(x = weather_median, y = -Inf, label = bin_label), vjust = -1.2, size = 3, color = "grey35", inherit.aes = FALSE) +
    annotate("text", x = Inf, y = hist_floor + 0.12 * y_span, label = "covariate distribution", hjust = 1.02, vjust = 0, size = 2.8, color = "grey45") +
    scale_fill_manual(values = c(low = "grey60", middle = "grey75", high = "grey60"), guide = "none") +
    labs(
      title = "Low-till duration effect by EDD, Jun-Sep bin",
      subtitle = "Pooled binning estimator with EDD bins and additive controls for the other weather variables.",
      x = paste0(weather_labels[[weather_id]], " (raw ", weather_units[[weather_id]], " units)"),
      y = y_label,
      caption = paste("Reference duration: clean conventional-till controls. Fixed effects: field and year. Weights: field size.", exclude_years_note)
    ) +
    coord_cartesian(ylim = c(hist_floor, y_ceiling + 0.08 * y_span), clip = "off") +
    theme_binning()
}

p <- make_single_weather_plot("EDD_6_9")
single_png <- file.path(plot_dir, paste0("corn_linear_duration_binning_edd69_with_other_weather_controls_cme_EDD_6_9_", scale_suffix, ".png"))
single_pdf <- file.path(plot_dir, paste0("corn_linear_duration_binning_edd69_with_other_weather_controls_cme_EDD_6_9_", scale_suffix, ".pdf"))
ggsave(single_png, p, width = 8.5, height = 5.2, dpi = 320)
ggsave(single_pdf, p, width = 8.5, height = 5.2)

plot_manifest <- data.table(
  plot_type = c("single_weather", "single_weather"),
  file = c(single_png, single_pdf)
)
fwrite(plot_manifest, file.path(plot_dir, paste0("corn_linear_duration_binning_edd69_with_other_weather_controls_cme_plot_manifest_", scale_suffix, ".csv")))

summary_sections <- c(
  "## Conditional Low-Till Duration Effects by Weather Bin",
  paste(capture.output(print(effects[, .(
    weather_var, weather_bin, weather_median, observations,
    estimate_log_points_per_year, std_error,
    pct_effect_per_year, pct_effect_low_per_year, pct_effect_high_per_year
  )])), collapse = "\n"),
  "",
  "## Additive Other-Weather Control Terms",
  paste(capture.output(print(additive_control_terms)), collapse = "\n"),
  "",
  "## Model Summaries",
  unlist(lapply(names(model_list), function(nm) {
    c(
      paste0("### ", nm),
      paste(capture.output(summary(model_list[[nm]])), collapse = "\n"),
      ""
    )
  }))
)

report_lines <- c(
  "# Binning Estimator for Conditional Low-Till Duration Effects on corn Yield (EDD only, With Other Weather Controls)",
  "",
  "## Design",
  "- Outcome: log corn yield in bushels per acre (`log_y`).",
  "- Treatment: continuous low-till duration in years (`low_till_duration`).",
  "- Control/reference duration: clean conventional-till observations with no prior low-till exposure (`till == 0` and `till_1_count == 0`).",
  "- Low-till observations: current low-till in an uninterrupted low-till spell (`till == 1` and `till_1_count == till_1_streak`).",
  paste0("- Excluded years: ", if (length(exclude_years) == 0) "none." else paste(exclude_years, collapse = ", "), "."),
  "- Focal weather moderator: `EDD_6_9` only.",
  "- Additive weather controls: `GDD_4_5`, `ppt_4_5`, `GDD_6_9`, and `ppt_6_9`.",
  "- Full weather sample: no focal-weather trimming is applied before constructing tercile bins and fitting the estimator.",
  "- Pooled estimator: the focal EDD moderator is split into low/middle/high tercile bins using the full focal-weather sample.",
  "- Model specification: control-group bin intercept dummies, bin-specific low-till duration effects, bin-specific focal-weather slopes, bin-specific duration-by-focal-weather interactions, and additive controls for the other weather variables.",
  "- Fixed effects: field and year.",
  paste0("- Standard errors clustered by `", cluster_var, "`."),
  "- Weights: field size.",
  "- Visualization note: the bottom histogram band shows only the middle 95% of the focal-weather distribution (2.5th-97.5th percentiles) to keep extreme tails from dominating the plot scale.",
  "",
  "## Interpretation",
  "- `estimate_log_points_per_year` is the conditional marginal effect of one additional year of continuous low-till duration at the median EDD value within that bin.",
  "- `pct_effect_per_year` converts that log-point estimate to percent using `100 * (exp(estimate) - 1)`.",
  "",
  "## Sample Counts",
  paste(capture.output(print(sample_counts)), collapse = "\n"),
  "",
  summary_sections,
  ""
)

writeLines(report_lines, file.path(output_dir, "corn_linear_duration_weather_binning_edd69_with_other_weather_controls_report.md"))

message("Done. Outputs written to ", output_dir)
message("Figure files written to ", plot_dir)
