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

root_dir <- find_root_dir()
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(
  root_dir,
  "output_linear_1999_2023",
  "1_corn_duration_bins_field_year_high_till_initially"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

crop <- "corn"
target_rows <- as.integer(Sys.getenv("TARGET_ROWS", unset = "80000000"))
sample_seed <- as.integer(Sys.getenv("SAMPLE_SEED", unset = "42"))

message("Root directory: ", root_dir)
message("Crop: ", crop)
message("Target rows: ", target_rows)

field_dir <- file.path(data_dir,
                       "Field_yield_tillage",
                       "processed_smooth_1999_2023")
pattern <- paste0("^", crop, "_yield_tillage_[A-Z]{2}_1999_2023\\.csv$")
field_files <- list.files(field_dir, pattern = pattern, full.names = TRUE)

if (length(field_files) == 0) {
  stop("No field files found for crop = ", crop)
}

field_cols <- c(
  "OBJECTID",
  "yield",
  "size",
  "year",
  "till",
  "GEOID",
  "till_1_count",
  "till_1_streak",
  "till_0_count",
  "till_0_streak"
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

dt <- dt[is.finite(yield_bu) &
           yield_bu > 0 & is.finite(size) & size > 0]

omit.fields = c(unique(dt[year == 1999 &
                            till_1_count == 1]$OBJECTID), unique(dt[year == 2000 &
                                                                      till_1_count == 2]$OBJECTID))

dt <- dt[!OBJECTID %in% omit.fields]

dt[, log_y := log(yield_bu)]

# Clean controls: never exposed to low-till up to the current year.
dt[, clean_control := till == 0L & till_1_count == 0L]

# Continuous adopters: current low-till and current streak equals cumulative count.
dt[, continuous_low_till := till == 1L &
     till_1_count == till_1_streak]

dt_duration <- dt[clean_control | continuous_low_till]

dt_duration[, till_duration_bin := fifelse(clean_control,
                                           "0",
                                           fifelse(
                                             till_1_streak <= 4L,
                                             "1_4",
                                             fifelse(
                                               till_1_streak <= 9L,
                                               "5_9",
                                               fifelse(
                                                 till_1_streak <= 14L,
                                                 "10_14",
                                                 fifelse(till_1_streak <= 19L, "15_19", "20+")
                                               )
                                             )
                                           ))]

dt_duration[, till_duration_bin := factor(till_duration_bin,
                                          levels = c("0", "1_4", "5_9", "10_14", "15_19", "20+"))]

climate_cols <- c("FIPS",
                  "year",
                  "GDD_4_5",
                  "ppt_4_5",
                  "GDD_6_9",
                  "ppt_6_9",
                  "EDD_6_9")
dt_climate <- fread(file.path(data_dir, "PRISM", "PRISM_climate_1999_2023.csv"),
                    select = climate_cols)
setnames(dt_climate, c("ppt_4_5", "ppt_6_9"), c("PPT_4_5", "PPT_6_9"))
dt_climate[, `:=`(
  FIPS = as.integer(FIPS),
  year = as.integer(year),
  GDD_4_5 = as.numeric(GDD_4_5),
  PPT_4_5 = as.numeric(PPT_4_5),
  GDD_6_9 = as.numeric(GDD_6_9),
  PPT_6_9 = as.numeric(PPT_6_9),
  EDD_6_9 = as.numeric(EDD_6_9)
)]

dt_duration <- merge(dt_duration,
                     dt_climate,
                     by = c("FIPS", "year"),
                     all.x = TRUE)
dt_duration <- dt_duration[is.finite(GDD_4_5) &
                             is.finite(PPT_4_5) &
                             is.finite(GDD_6_9) &
                             is.finite(PPT_6_9) &
                             is.finite(EDD_6_9)]

sample_counties_by_state <- function(x, n_rows, seed) {
  if (is.na(n_rows) || n_rows <= 0 || n_rows >= nrow(x)) {
    return(copy(x))
  }
  
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
    if (length(keep_fips) == 0) {
      keep_fips <- county_rows[1, FIPS]
    }
    sampled_parts[[i]] <- x[FIPS %in% keep_fips]
  }
  
  rbindlist(sampled_parts, use.names = TRUE)
}

dt_sample <- sample_counties_by_state(dt_duration, target_rows, sample_seed)
dt_sample[, state := as.factor((FIPS - (FIPS %% 1000L)) / 1000L)]

dt_irrigation <- fread(file.path(
  data_dir,
  "irrigation",
  "field_irrigation_ever_only_irrigated.csv"
))
dt_sample <- dt_sample[!dt_sample$OBJECTID %in% dt_irrigation$OBJECTID, ]

covars <- "GDD_4_5 + PPT_4_5 + GDD_6_9 + PPT_6_9 + EDD_6_9"

model_logy_tillage_duration_field_year_climate <- feols(
  as.formula(
    paste0(
      "log_y ~ i(till_duration_bin, ref = '0') + ",
      covars,
      " | OBJECTID + year"
    )
  ),
  cluster = ~ FIPS,
  weights = ~ size,
  lean = TRUE,
  mem.clean = TRUE,
  data = dt_sample
)

model_summary <- summary(model_logy_tillage_duration_field_year_climate)
print(model_summary)

term_table <- as.data.table(model_summary$coeftable, keep.rownames = "term")
setnames(
  term_table,
  old = c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
  new = c("estimate", "std_error", "t_value", "p_value"),
  skip_absent = TRUE
)
term_table[, model := "field_year_climate"]
setcolorder(term_table, c("model", "term"))

sample_counts <- dt_sample[, .(
  observations = .N,
  fields = uniqueN(OBJECTID),
  counties = uniqueN(FIPS),
  acres = sum(size, na.rm = TRUE)
), by = till_duration_bin]
setorder(sample_counts, till_duration_bin)

fwrite(
  sample_counts,
  file.path(
    output_dir,
    "corn_fe_duration_bins_clean_climate_sample_counts.csv"
  )
)
fwrite(
  term_table,
  file.path(output_dir, "corn_fe_duration_bins_clean_climate_terms.csv")
)
writeLines(
  capture.output(model_summary),
  file.path(
    output_dir,
    "corn_fe_duration_bins_clean_climate_summary.txt"
  )
)
saveRDS(
  list(field_year_climate = model_logy_tillage_duration_field_year_climate),
  file.path(
    output_dir,
    "corn_fe_duration_bins_clean_climate_models.rds"
  )
)

message("Done. Outputs written to ", output_dir)
