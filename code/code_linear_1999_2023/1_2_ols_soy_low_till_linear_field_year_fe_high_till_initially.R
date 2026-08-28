## =============================== PORTABLE PATHS ==============================
## Set ROOT_OVERRIDE to the folder that contains "data/" and leave everything else
## alone. If you leave it as "", the script locates the project automatically.
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

extract_terms <- function(model) {
  out <- as.data.table(coeftable(model), keep.rownames = "term")
  setnames(
    out,
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
    c("estimate", "std_error", "t_value", "p_value")
  )
  out[, `:=`(
    pct_effect = (exp(estimate) - 1) * 100,
    pct_effect_low = (exp(estimate - 1.96 * std_error) - 1) * 100,
    pct_effect_high = (exp(estimate + 1.96 * std_error) - 1) * 100
  )]
  out[]
}

root_dir <- find_root_dir()
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "output_linear_1999_2023", "1_soy_linear_field_year_high_till_initially")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

field_dir <- file.path(data_dir, "Field_yield_tillage", "processed_smooth_1999_2023")
field_files <- list.files(
  field_dir,
  pattern = "^soy_yield_tillage_[A-Z]{2}_1999_2023\\.csv$",
  full.names = TRUE
)

if (length(field_files) == 0) {
  stop("No soybean field-level tillage/yield files found in ", field_dir)
}

field_cols <- c(
  "OBJECTID", "yield", "size", "year", "till", "GEOID",
  "till_1_count", "till_1_streak", "till_0_count", "till_0_streak"
)

message("Loading soybean field files: ", length(field_files))
dt <- rbindlist(
  lapply(field_files, fread, select = field_cols),
  use.names = TRUE
)

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

dt <- dt[
  is.finite(yield_bu) & yield_bu > 0 &
    is.finite(size) & size > 0 &
    !is.na(OBJECTID) & !is.na(FIPS) & !is.na(year) & !is.na(till)
]

omit.fields = c(unique(dt[year == 1999 &
                            till_1_count == 1]$OBJECTID), unique(dt[year == 2000 &
                                                                      till_1_count == 2]$OBJECTID))

dt <- dt[!OBJECTID %in% omit.fields]

dt[, log_yield := log(yield_bu)]

# Clean controls: conventional observations with no previous low-till exposure.
dt[, clean_control := till == 0L & till_1_count == 0L]

# Continuous low-till: current low-till observations in an uninterrupted low-till spell.
dt[, continuous_low_till := till == 1L & till_1_count == till_1_streak]

sample_before_eligibility_filter <- dt[, .(
  sample = "full_cleaned_panel",
  observations = .N,
  fields = uniqueN(OBJECTID),
  counties = uniqueN(FIPS),
  years = uniqueN(year)
)]

dt <- dt[clean_control | continuous_low_till]

sample_before_weather <- dt[, .(
  sample = "clean_control_or_continuous_low_till",
  observations = .N,
  fields = uniqueN(OBJECTID),
  counties = uniqueN(FIPS),
  years = uniqueN(year)
)]

# Treatment: current duration of continuous low-till exposure.
# Clean controls receive duration 0; continuous low-till observations use the current low-till streak.
dt[, low_till_duration := fifelse(clean_control, 0L, till_1_streak)]
dt <- dt[is.finite(low_till_duration) & low_till_duration >= 0]

weather_cols <- c("FIPS", "year", "GDD_4_5", "ppt_4_5", "GDD_6_9", "ppt_6_9", "EDD_6_9")
weather <- fread(
  file.path(data_dir, "PRISM", "PRISM_climate_1999_2023.csv"),
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

dt_irrigation <- fread(file.path(data_dir, "irrigation", "field_irrigation_ever_only_irrigated.csv"))
dt <- dt[!dt$OBJECTID %in% dt_irrigation$OBJECTID, ]

sample_counts <- rbindlist(list(
  sample_before_eligibility_filter,
  sample_before_weather,
  dt[, .(
    sample = "estimation_sample_after_weather_join",
    observations = .N,
    fields = uniqueN(OBJECTID),
    counties = uniqueN(FIPS),
    years = uniqueN(year)
  )]
), use.names = TRUE)

eligibility_counts <- dt[, .(
  observations = .N,
  fields = uniqueN(OBJECTID)
), by = .(eligibility_group = fifelse(clean_control, "clean_control", "continuous_low_till"))]

message("Estimation observations: ", format(nrow(dt), big.mark = ","))
message("Estimation fields: ", format(uniqueN(dt$OBJECTID), big.mark = ","))

ols_field_year <- feols(
  log_yield ~ low_till_duration + GDD_4_5 + PPT_4_5 + GDD_6_9 + PPT_6_9 + EDD_6_9 |
    OBJECTID + year,
  data = dt,
  weights = ~size,
  cluster = ~FIPS,
  lean = TRUE,
  mem.clean = TRUE
)

print(summary(ols_field_year))

term_table <- extract_terms(ols_field_year)

fwrite(sample_counts, file.path(output_dir, "soy_sample_counts.csv"))
fwrite(eligibility_counts, file.path(output_dir, "soy_eligibility_counts.csv"))
fwrite(term_table, file.path(output_dir, "soy_ols_field_year_terms.csv"))
saveRDS(ols_field_year, file.path(output_dir, "soy_ols_field_year_model.rds"))

message("Outputs written to: ", output_dir)
