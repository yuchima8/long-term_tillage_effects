library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[[1]])
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }
  getwd()
}

find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  for (i in seq_len(10)) {
    if (
      dir.exists(file.path(cur, "data", "Field_yield_tillage", "processed_smooth_1999_2023")) &&
      file.exists(file.path(cur, "data", "PRISM", "PRISM_climate_1999_2023.csv"))
    ) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Could not locate project root from ", start_dir)
}

make_weather_histograms <- function(
    root_dir,
    crop,
    weather_order_levels,
    weather_labels,
    n_bins = 35L,
    trim_probs = c(0.025, 0.975)
) {
  field_dir <- file.path(root_dir, "data", "Field_yield_tillage", "processed_smooth_1999_2023")
  field_files <- list.files(
    field_dir,
    pattern = paste0("^", crop, "_yield_tillage_[A-Z]{2}_1999_2023\\.csv$"),
    full.names = TRUE
  )
  if (length(field_files) == 0) stop("No field files found for crop = ", crop)
  
  field_cols <- c(
    "OBJECTID", "yield", "size", "year", "till", "GEOID",
    "till_1_count", "till_1_streak"
  )
  
  message("Building weather histograms from field files: ", length(field_files))
  dt <- rbindlist(lapply(field_files, fread, select = field_cols), use.names = TRUE)
  setnames(dt, "GEOID", "FIPS")
  dt[, `:=`(
    OBJECTID = as.character(OBJECTID),
    FIPS = as.integer(FIPS),
    year = as.integer(year),
    till = as.integer(till),
    till_1_count = as.integer(till_1_count),
    till_1_streak = as.integer(till_1_streak),
    size = as.numeric(size),
    yield_bu = as.numeric(yield) / 1000
  )]
  
  dt <- dt[
    is.finite(yield_bu) & yield_bu > 0 &
      is.finite(size) & size > 0 &
      !is.na(OBJECTID) & !is.na(FIPS) & !is.na(year) & !is.na(till)
  ]
  dt[, clean_control := till == 0L & till_1_count == 0L]
  dt[, continuous_low_till := till == 1L & till_1_count == till_1_streak]
  dt <- dt[clean_control | continuous_low_till]
  dt[, low_till_duration := fifelse(clean_control, 0L, till_1_streak)]
  dt <- dt[is.finite(low_till_duration) & low_till_duration >= 0]
  
  weather_cols <- c("FIPS", "year", "GDD_4_5", "ppt_4_5", "GDD_6_9", "ppt_6_9", "EDD_6_9")
  weather <- fread(
    file.path(root_dir, "data", "PRISM", "PRISM_climate_1999_2023.csv"),
    select = weather_cols
  )
  setnames(weather, c("ppt_4_5", "ppt_6_9"), c("PPT_4_5", "PPT_6_9"))
  weather[, `:=`(
    FIPS = as.integer(FIPS),
    year = as.integer(year),
    GDD_4_5 = as.numeric(GDD_4_5),
    PPT_4_5 = as.numeric(PPT_4_5),
    GDD_6_9 = as.numeric(GDD_6_9),
    PPT_6_9 = as.numeric(PPT_6_9),
    EDD_6_9 = as.numeric(EDD_6_9)
  )]
  
  dt <- merge(dt, weather, by = c("FIPS", "year"), all.x = TRUE, sort = FALSE)
  dt <- dt[
    is.finite(GDD_4_5) &
      is.finite(PPT_4_5) &
      is.finite(GDD_6_9) &
      is.finite(PPT_6_9) &
      is.finite(EDD_6_9)
  ]

  irrigation_path <- file.path(root_dir, "data", "irrigation", "field_irrigation_ever_only_irrigated.csv")
  if (file.exists(irrigation_path)) {
    irrigated_fields <- fread(irrigation_path, select = "OBJECTID")
    irrigated_fields[, OBJECTID := as.character(OBJECTID)]
    dt <- dt[!OBJECTID %in% irrigated_fields$OBJECTID]
  }
  
  out <- vector("list", length(weather_order_levels))
  for (i in seq_along(weather_order_levels)) {
    weather_var <- weather_order_levels[[i]]
    
    x <- dt[[weather_var]]
    x <- x[is.finite(x)]
    
    bounds <- stats::quantile(
      x,
      probs = trim_probs,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    
    if (!all(is.finite(bounds)) || bounds[[1]] >= bounds[[2]]) {
      stop("Could not compute middle 95% histogram bounds for ", weather_var)
    }
    
    x <- x[x >= bounds[[1]] & x <= bounds[[2]]]
    breaks <- seq(bounds[[1]], bounds[[2]], length.out = n_bins + 1L)
    hist_counts <- hist(x, breaks = breaks, plot = FALSE, include.lowest = TRUE, right = TRUE)$counts
    
    out[[i]] <- data.frame(
      weather_var = weather_var,
      weather_label = unname(weather_labels[[weather_var]]),
      xmin = head(breaks, -1L),
      xmax = tail(breaks, -1L),
      count = hist_counts
    )
  }
  
  bind_rows(out)
}

## ---------------------------------------------------------------------------
## WHERE ARE THE RESULT CSVs?
## Leave RESULT_DIR as "" and the script searches the usual places. Set it if the
## search fails -- typically because the code and the outputs live on different
## drives. Point it at the folder that CONTAINS the four result CSVs, e.g.
##   RESULT_DIR <- "/Volumes/T9/yuchi/Stanford/tillage_effect/output_linear_1999_2023/2_1_soy_low_till_duration_weather_binning_25_50_25"
## Equivalent alternatives, no editing required:
##   Rscript summarize_soybean_binning_estimator.R <that folder>
##   BINNING_RESULT_DIR=<that folder> Rscript summarize_soybean_binning_estimator.R
## ---------------------------------------------------------------------------
RESULT_DIR <- ""

args <- commandArgs(trailingOnly = TRUE)
if (!nzchar(RESULT_DIR)) RESULT_DIR <- Sys.getenv("BINNING_RESULT_DIR", unset = "")
base_dir <- if (nzchar(RESULT_DIR)) {
  normalizePath(RESULT_DIR, mustWork = TRUE)
} else if (length(args) >= 1) {
  normalizePath(args[[1]], mustWork = TRUE)
} else {
  script_dir()
}

## The project root is needed ONLY to build the weather histogram for the point
## plot, and that histogram is cached to CSV after the first run. So a missing
## data/ folder must not be fatal here -- it is checked later, at the point of use.
root_dir <- tryCatch(find_project_root(base_dir), error = function(e) NA_character_)

result_tag <- "25_50_25_with_other_weather_controls"
target_crop <- "soy"
target_result_dir <- "2_1_soy_low_till_duration_weather_binning_25_50_25"

resolve_result_files <- function(base_dir, root_dir) {
  prefix_i <- paste0(target_crop, "_linear_duration_weather_binning_", result_tag)
  suffixes <- c(effects = "_effects.csv",
                duration_slopes = "_duration_weather_slopes.csv",
                control_slopes = "_control_weather_slopes.csv",
                bin_stats = "_bin_stats.csv")

  ## Search order: the explicit/default base_dir, the script's own folder, then the
  ## conventional output folder relative to the project root AND relative to the
  ## code folder. The last of these is what makes the script work when data/ is
  ## absent, because it does not depend on find_project_root() having succeeded.
  candidate_dirs <- c(base_dir, script_dir())
  if (!is.na(root_dir)) {
    candidate_dirs <- c(candidate_dirs,
      file.path(root_dir, "output_linear_1999_2023", target_result_dir),
      file.path(root_dir, "output_linear_1999_2023",
                paste0(target_crop, "_low_till_duration_weather_binning_", result_tag)),
      file.path(root_dir, "output_linear_1999_2023",
                paste0("2_", target_crop, "_low_till_duration_weather_binning_", result_tag)))
  }
  for (start in unique(c(base_dir, script_dir()))) {
    up1 <- dirname(start)          # e.g. .../code
    up2 <- dirname(up1)            # e.g. the project root
    candidate_dirs <- c(candidate_dirs,
      file.path(up1, "output_linear_1999_2023", target_result_dir),
      file.path(up2, "output_linear_1999_2023", target_result_dir))
  }
  candidate_dirs <- unique(candidate_dirs)

  searched <- character()
  for (dir_i in candidate_dirs) {
    paths_i <- stats::setNames(file.path(dir_i, paste0(prefix_i, suffixes)), names(suffixes))
    if (!dir.exists(dir_i)) {
      searched <- c(searched, sprintf("  no such folder                  %s", dir_i))
      next
    }
    hit <- file.exists(paths_i)
    if (all(hit)) {
      message("Reading result CSVs from: ", normalizePath(dir_i))
      return(list(base_dir = normalizePath(dir_i), crop = target_crop,
                  prefix = prefix_i, paths = paths_i))
    }
    searched <- c(searched, sprintf("  folder exists, %d of 4 CSVs     %s", sum(hit), dir_i))
  }

  ## R truncates stop() messages at getOption("warning.length") = 1000 chars, which
  ## cut off the FIX instructions. message() is not truncated, so print the detail
  ## first and keep the stop() condition short.
  message(paste0(
       "Could not find the four result CSVs for crop = ", target_crop, ".\n\n",
       "Looking for these file names:\n  ",
       paste(paste0(prefix_i, suffixes), collapse = "\n  "), "\n\n",
       "Searched (in order):\n", paste(searched, collapse = "\n"), "\n\n",
       "FIX. Set RESULT_DIR near the top of this script to the folder holding those\n",
       "files -- normally  output_linear_1999_2023/", target_result_dir, "  -- or pass\n",
       "that folder as the first argument to Rscript.\n\n",
       "If none of the searched folders contains any of the CSVs, the upstream model\n",
       "has not been run (or its output was not copied to this machine). Run\n",
       "  2_1_ols_soy_low_till_duration_weather_binning_25_50_25.R\n",
       "first; this script only summarises that script's output."))
  stop("Result CSVs not found for crop = ", target_crop,
       " -- see the diagnostic printed above. Set RESULT_DIR at the top of this script.")
}

resolved <- resolve_result_files(base_dir, root_dir)
base_dir <- resolved$base_dir
crop <- resolved$crop
crop_label <- ifelse(crop %in% c("soy", "soybean"), "soybean", "corn")
prefix <- resolved$prefix

effects_path <- resolved$paths[["effects"]]
duration_slopes_path <- resolved$paths[["duration_slopes"]]
control_slopes_path <- resolved$paths[["control_slopes"]]
bin_stats_path <- resolved$paths[["bin_stats"]]

out_dir <- base_dir

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

effects <- read.csv(effects_path, check.names = FALSE)
duration_slopes <- read.csv(duration_slopes_path, check.names = FALSE)
control_slopes <- read.csv(control_slopes_path, check.names = FALSE)
bin_stats <- read.csv(bin_stats_path, check.names = FALSE)

for (optional_col in c("trimmed_sample_p05", "trimmed_sample_p95")) {
  if (!optional_col %in% names(bin_stats)) {
    bin_stats[[optional_col]] <- NA_real_
  }
}

weather_labels <- c(
  GDD_4_5 = "GDD, Apr-May",
  PPT_4_5 = "PPT, Apr-May",
  GDD_6_9 = "GDD, Jun-Sep",
  PPT_6_9 = "PPT, Jun-Sep",
  EDD_6_9 = "EDD, Jun-Sep"
)

weather_units <- c(
  GDD_4_5 = "degree days",
  PPT_4_5 = "mm",
  GDD_6_9 = "degree days",
  PPT_6_9 = "mm",
  EDD_6_9 = "degree days"
)

weather_order_levels <- c("GDD_4_5", "PPT_4_5", "GDD_6_9", "PPT_6_9", "EDD_6_9")
weather_bin_levels <- c("low", "middle", "high")
weather_panel_labels <- paste0(
  unname(weather_labels[weather_order_levels]),
  "\n(",
  unname(weather_units[weather_order_levels]),
  ")"
)

fmt_num <- function(x, digits = 1) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_effect_ci <- function(effect, low, high) {
  paste0(fmt_num(effect, 3), " (", fmt_num(low, 3), ", ", fmt_num(high, 3), ")")
}

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ ""
  )
}

summary_long <- effects %>%
  mutate(
    weather_var = factor(weather_var, levels = weather_order_levels),
    weather_label = unname(weather_labels[as.character(weather_var)]),
    weather_unit = unname(weather_units[as.character(weather_var)]),
    weather_bin = factor(weather_bin, levels = weather_bin_levels),
    bin_label = recode(as.character(weather_bin),
                       low = "Low",
                       middle = "Middle",
                       high = "High"
    ),
    bin_range = paste0(
      fmt_num(weather_min, 1), "-",
      fmt_num(weather_max, 1), " ",
      weather_unit
    ),
    median_with_unit = paste0(fmt_num(weather_median, 1), " ", weather_unit),
    t_value = estimate_log_points_per_year / std_error,
    p_value = 2 * pnorm(abs(t_value), lower.tail = FALSE),
    significant_05 = p_value < 0.05,
    significance = sig_stars(p_value),
    effect_ci = fmt_effect_ci(pct_effect_per_year, pct_effect_low_per_year, pct_effect_high_per_year)
  ) %>%
  arrange(weather_var, weather_bin)

summary_table <- summary_long %>%
  transmute(
    `Weather variable` = weather_label,
    `Weather bin` = bin_label,
    `Median weather value` = median_with_unit,
    `Bin range` = bin_range,
    `Observations` = observations,
    `Fields` = fields,
    `Counties` = counties,
    `Control observations` = control_observations,
    `Low-till observations` = low_till_observations,
    `Mean low-till duration` = round(mean_low_till_duration, 2),
    `Median low-till duration` = round(median_low_till_duration, 2),
    `Low-till duration effect per year (%; 95% CI)` = effect_ci,
    `p-value` = signif(p_value, 3),
    `Significant at 5%` = significant_05,
    `Significance` = significance
  )

write.csv(
  summary_table,
  file.path(out_dir, paste0(crop, "_binning_estimator_duration_effect_summary_table.csv")),
  row.names = FALSE
)

writeLines(
  c(
    "| Weather variable | Weather bin | Median weather value | Bin range | Observations | Fields | Counties | Control observations | Low-till observations | Mean low-till duration | Median low-till duration | Low-till duration effect per year (%; 95% CI) | p-value | Significant at 5% | Significance |",
    "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    apply(summary_table, 1, function(row) {
      paste0("| ", paste(row, collapse = " | "), " |")
    })
  ),
  con = file.path(out_dir, paste0(crop, "_binning_estimator_duration_effect_summary_table.md"))
)

main_numbers <- summary_long %>%
  group_by(weather_label) %>%
  slice_max(order_by = pct_effect_per_year, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    `Weather variable` = weather_label,
    `Largest effect bin` = bin_label,
    `Duration effect per year (%)` = round(pct_effect_per_year, 3),
    `95% CI low (%)` = round(pct_effect_low_per_year, 3),
    `95% CI high (%)` = round(pct_effect_high_per_year, 3),
    `p-value` = signif(p_value, 3),
    `Significance` = significance
  )

write.csv(
  main_numbers,
  file.path(out_dir, paste0(crop, "_binning_estimator_main_duration_numbers.csv")),
  row.names = FALSE
)

slope_summary <- duration_slopes %>%
  select(
    weather_var, weather_bin, bin_order,
    slope_log_points_per_year_per_weather_unit,
    slope_std_error
  ) %>%
  left_join(
    control_slopes %>%
      select(weather_var, weather_bin, control_weather_slope, control_weather_slope_std_error),
    by = c("weather_var", "weather_bin")
  ) %>%
  mutate(
    weather_var = factor(weather_var, levels = weather_order_levels),
    weather_label = unname(weather_labels[as.character(weather_var)]),
    weather_unit = unname(weather_units[as.character(weather_var)]),
    weather_bin = factor(weather_bin, levels = weather_bin_levels),
    bin_label = recode(as.character(weather_bin),
                       low = "Low",
                       middle = "Middle",
                       high = "High"
    ),
    duration_weather_slope_pct = 100 * slope_log_points_per_year_per_weather_unit,
    duration_weather_slope_pct_se = 100 * slope_std_error,
    duration_weather_slope_p_value = 2 * pnorm(
      abs(slope_log_points_per_year_per_weather_unit / slope_std_error),
      lower.tail = FALSE
    ),
    control_weather_slope_pct = 100 * control_weather_slope,
    control_weather_slope_pct_se = 100 * control_weather_slope_std_error,
    control_weather_slope_p_value = 2 * pnorm(
      abs(control_weather_slope / control_weather_slope_std_error),
      lower.tail = FALSE
    )
  ) %>%
  arrange(weather_var, weather_bin) %>%
  transmute(
    `Weather variable` = weather_label,
    `Weather bin` = bin_label,
    `Weather unit` = weather_unit,
    `Control weather slope (% yield per weather unit)` = round(control_weather_slope_pct, 5),
    `Control weather slope SE` = round(control_weather_slope_pct_se, 5),
    `Control weather slope p-value` = signif(control_weather_slope_p_value, 3),
    `Duration-weather slope (% yield per year per weather unit)` = round(duration_weather_slope_pct, 5),
    `Duration-weather slope SE` = round(duration_weather_slope_pct_se, 5),
    `Duration-weather slope p-value` = signif(duration_weather_slope_p_value, 3),
    `Duration-weather slope significance` = sig_stars(duration_weather_slope_p_value)
  )

write.csv(
  slope_summary,
  file.path(out_dir, paste0(crop, "_binning_estimator_weather_slope_summary.csv")),
  row.names = FALSE
)

bin_summary <- bin_stats %>%
  mutate(
    weather_var = factor(weather_var, levels = weather_order_levels),
    weather_label = unname(weather_labels[as.character(weather_var)]),
    weather_bin = factor(weather_bin, levels = weather_bin_levels),
    bin_label = recode(as.character(weather_bin),
                       low = "Low",
                       middle = "Middle",
                       high = "High"
    )
  ) %>%
  arrange(weather_var, weather_bin) %>%
  transmute(
    `Weather variable` = weather_label,
    `Weather bin` = bin_label,
    `Observations` = observations,
    `Fields` = fields,
    `Counties` = counties,
    `Control observations` = control_observations,
    `Low-till observations` = low_till_observations,
    `Mean low-till duration` = round(mean_low_till_duration, 2),
    `Median low-till duration` = round(median_low_till_duration, 2),
    `Weather min` = round(weather_min, 2),
    `Weather median` = round(weather_median, 2),
    `Weather max` = round(weather_max, 2),
    `Sample p05` = round(trimmed_sample_p05, 2),
    `Sample p95` = round(trimmed_sample_p95, 2)
  )

write.csv(
  bin_summary,
  file.path(out_dir, paste0(crop, "_binning_estimator_bin_sample_summary.csv")),
  row.names = FALSE
)

point_data <- summary_long %>%
  mutate(
    weather_panel = factor(
      paste0(weather_label, "\n(", weather_unit, ")"),
      levels = weather_panel_labels
    ),
    value_label = paste0(sprintf("%.2f", pct_effect_per_year), "%"),
    label_x = weather_median + 0.035 * (weather_max - weather_min),
    label_y = pct_effect_per_year
  )

y_range <- range(c(point_data$pct_effect_low_per_year, point_data$pct_effect_high_per_year), na.rm = TRUE)
y_pad <- diff(y_range) * 0.12
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.1
distribution_y <- y_range[[1]] - y_pad
histogram_height <- y_pad * 0.55

histogram_cache_path <- file.path(out_dir, paste0(crop, "_weather_covariate_histograms_middle95.csv"))
if (file.exists(histogram_cache_path)) {
  histogram_data_raw <- read.csv(histogram_cache_path, check.names = FALSE)
} else {
  if (is.na(root_dir)) {
    stop("The weather histogram for the point plot has not been cached yet, and it is\n",
         "built from the raw field panel, which was not found. Either:\n",
         "  - make data/Field_yield_tillage/processed_smooth_1999_2023/ and\n",
         "    data/PRISM/PRISM_climate_1999_2023.csv reachable from this folder, or\n",
         "  - copy an existing ", basename(histogram_cache_path), " into\n",
         "    ", out_dir)
  }
  histogram_data_raw <- make_weather_histograms(
    root_dir = root_dir,
    crop = crop,
    weather_order_levels = weather_order_levels,
    weather_labels = weather_labels,
    n_bins = 35L,
    trim_probs = c(0.025, 0.975)
  )
  write.csv(histogram_data_raw, histogram_cache_path, row.names = FALSE)
}

histogram_data <- histogram_data_raw %>%
  mutate(
    weather_panel = factor(
      paste0(weather_label, "\n(", weather_units[weather_var], ")"),
      levels = weather_panel_labels
    )
  ) %>%
  group_by(weather_panel) %>%
  mutate(hist_scaled = count / max(count, na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(
    weather_panel,
    xmin,
    xmax,
    ymin = distribution_y,
    ymax = distribution_y + histogram_height * hist_scaled
  )

boundary_data <- point_data %>%
  filter(as.character(weather_bin) %in% c("low", "middle")) %>%
  transmute(
    weather_panel,
    boundary = weather_max
  )

p_points <- ggplot(point_data, aes(x = weather_median, y = pct_effect_per_year)) +
  geom_hline(yintercept = 0, color = "gray50", linewidth = 0.4) +
  geom_vline(
    data = boundary_data,
    aes(xintercept = boundary),
    linetype = "22",
    color = "gray40",
    linewidth = 0.45
  ) +
  geom_rect(
    data = histogram_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "gray80",
    color = "gray55",
    linewidth = 0.25
  ) +
  geom_linerange(aes(ymin = pct_effect_low_per_year, ymax = pct_effect_high_per_year),
                 linewidth = 0.55,
                 color = "#3D6C91"
  ) +
  geom_point(shape = 21, size = 2.9, color = "gray15", fill = "#5B8DB8", stroke = 0.35) +
  geom_label(aes(x = label_x, y = label_y, label = value_label),
             hjust = 0,
             size = 2.75,
             fontface = "bold",
             linewidth = 0,
             label.padding = unit(0.08, "lines"),
             fill = scales::alpha("white", 0.82),
             color = "gray15"
  ) +
  facet_wrap(~ weather_panel, ncol = 5, scales = "free_x") +
  coord_cartesian(ylim = c(distribution_y - y_pad * 0.15, y_range[[2]] + y_pad), clip = "off") +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Yield effect per low-till year (%)",
    caption = paste(
      "Points and bars: duration effect per low-till year with 95% CI, plotted at the bin median.",
      "Dashed lines: bin boundaries.",
      "\nGray band: distribution of the weather variable (middle 95% of observations), rescaled to fit;",
      "it is NOT on the y-axis scale."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.25),
    panel.border = element_rect(color = "gray35", fill = NA, linewidth = 0.55),
    strip.text = element_text(face = "bold", size = 12, margin = margin(b = 5)),
    strip.background = element_rect(fill = "gray94", color = "gray35", linewidth = 0.55),
    axis.text.x = element_text(color = "gray15"),
    axis.text.y = element_text(color = "gray15"),
    plot.caption = element_text(color = "gray35", size = 7.6, hjust = 0),
    plot.margin = margin(10, 12, 8, 18)
  )

ggsave(file.path(out_dir, paste0(crop, "_binning_estimator_duration_effect_by_bin.png")), p_points, width = 12, height = 3.6, dpi = 350)
ggsave(file.path(out_dir, paste0(crop, "_binning_estimator_duration_effect_by_bin.pdf")), p_points, width = 11.2, height = 4.7, device = "pdf")

message("Summary files written to: ", out_dir)
