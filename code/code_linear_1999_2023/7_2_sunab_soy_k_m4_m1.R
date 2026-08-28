## =============================== PORTABLE PATHS ==============================
## Set ROOT_OVERRIDE to the folder that contains "data/" and leave everything else
## alone. If you leave it as "", the script locates the project automatically.
##
##   Google Drive (Mac)      ROOT_OVERRIDE <- "/Users/<you>/Library/CloudStorage/GoogleDrive-<you>@stanford.edu/My Drive/Research/Tillage_trend/Codex"
##   Google Drive (Windows)  ROOT_OVERRIDE <- "G:/My Drive/Research/Tillage_trend/Codex"
##   Cluster / anywhere      ROOT_OVERRIDE <- "/scratch/users/<you>/Tillage_effects"
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

###############################################################################
##  Sun-Abraham interaction-weighted (IW) event study -- SOYBEAN
##  Reference period: k = -4:-1
##
##  Robustness check for the two-way fixed-effects duration models. With staggered
##  adoption and effects that grow with exposure, the TWFE coefficient can place
##  negative weight on some treated field-years. The IW estimator saturates on
##  cohort x relative-period, so an already-treated cohort is never used as a
##  control for a later-treated one, and aggregates with cohort area shares, which
##  are non-negative by construction.
##
##  For an absorbing treatment, event time and duration are the same quantity:
##      adoption year a = min(year - till_1_streak + 1)
##      k = year - a    =>    duration = k + 1
##  so the relative-period profile IS a duration profile.
##
##  ---------------------------------------------------------------------------
##  HOW TO RUN
##  ---------------------------------------------------------------------------
##  1. Set ROOT below to the folder that contains "data/".
##       Google Drive desktop (Mac):  "/Users/<you>/Google Drive/Tillage_effects"
##       Google Drive desktop (Win):  "G:/My Drive/Tillage_effects"
##       Colab: run  library(googledrive)  or mount the drive, then point ROOT at it.
##  2. Expected layout under ROOT:
##       data/Field_yield_tillage/processed_smooth_1999_2023/soy_yield_tillage_XX_1999_2023.csv
##       data/PRISM/PRISM_climate_1999_2023.csv
##       data/irrigation/field_irrigation_ever_only_irrigated.csv
##  3. MEMORY. Measured peak on the full nine-state run: 20 GB (corn), 16 GB (soy),
##     about 2.5 minutes on 8 cores. Plan for ~24 GB. If you have
##     less, set STATES to a few states and/or FIELD_FRAC below 1. A quick smoke
##     test that runs in a couple of minutes on a laptop:
##       STATES <- c("IA");  FIELD_FRAC <- 0.25
###############################################################################

##  ---------------------------------------------------------------------------
##  THIS VARIANT: reference k = -4:-1
##  ---------------------------------------------------------------------------
##  Identical to 7_sunab_soy.R except that the omitted category is the four
##  years immediately before adoption instead of the whole pre-period.
##
##  What this buys you: the distant pre-adoption cells (<=-11, -10:-8, -7:-5) are
##  no longer absorbed into the reference, so they are estimated and you can read
##  them as a pre-trend test. Under the default k<0 reference every pre-period
##  cell is the reference, so no pre-trend test is possible.
##
##  EXPECTED: the rows -4:-3, -2 and -1 come back as exactly 0 with a blank CI and
##  is_reference = TRUE. That is not a bug and not a finding -- those three bins ARE
##  the omitted category, so they have no coefficient to report. Only <=-11,
##  -10:-8 and -7:-5 are informative about pre-trends here.
##
##  CAUTION on interpretation: k = -1 sits materially above a {-4:-3} baseline in
##  both crops (+0.49% corn, +1.53% soy, both significant -- an Ashenfelter-style
##  dip in reverse). Folding it into the reference shifts the whole post-period
##  profile up relative to the k<0 run. Compare levels across references only with
##  that in mind; use 7_sunab_soy.R (k<0) for the TWFE comparison, and
##  REFERENCE <- "k=-4:-3" if you want to see the k = -1 spike itself.
##  ---------------------------------------------------------------------------

## ============================ USER SETTINGS =================================
ROOT       <- ""          # leave "" to auto-detect; see ROOT_OVERRIDE above
OUTDIR     <- file.path(ROOT, "7_output_sunab")
CROP       <- "soy"

## Reference period: which relative-period cells go into the omitted category.
##   "k<0"      all pre-adoption cells   -- reproduces the TWFE reference group
##                                          (every field-year at zero duration);
##                                          use this to compare against TWFE.
##   "k=-4:-1"  four-year window         -- leaves the distant pre-cells estimable
##   "k=-4:-3"  two-year window          -- the only option that leaves k=-1
##                                          estimable, so it can show the k=-1 spike
##   "k=-2:-1"  immediate two years      -- NOT recommended: k = -1 is anomalously
##                                          high in both crops and contaminates it
REFERENCE  <- "k=-4:-1"

FE         <- "field_year"   # "field_year" or "field_year_state"
STATES     <- NULL           # NULL = all nine; else e.g. c("IA","IL")
FIELD_FRAC <- 1              # 1 = all fields; else a fraction, e.g. 0.25
SEED       <- 42
MIN_CELL   <- 5000           # cohort x period cells smaller than this join the reference
## ============================================================================

pkgs <- c("data.table", "fixest")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing package(s): ", paste(miss, collapse = ", "),
                       "\nInstall with: install.packages(c(",
                       paste(sprintf('"%s"', miss), collapse = ", "), "))")
library(data.table); library(fixest)
setDTthreads(0L)                     # 0 = use all available cores
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

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

if (!nzchar(ROOT)) ROOT <- find_root_dir()
OUTDIR <- file.path(ROOT, "output_sunab")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
DATA <- file.path(ROOT, "data")
stopifnot(dir.exists(DATA))
ST <- c("17"="IL","18"="IN","19"="IA","26"="MI","27"="MN",
        "29"="MO","39"="OH","46"="SD","55"="WI")

## relative-period bins. -2 is kept separate from -1 so that the k = -1 cell can be
## estimated under the k=-4:-3 reference; do NOT merge them (they differ materially).
RL <- c("<=-11","-10:-8","-7:-5","-4:-3","-2","-1","0","1:3","4:8","9:13","14:18",">=19")
REFSETS <- list("k<0"     = c("<=-11","-10:-8","-7:-5","-4:-3","-2","-1"),
                "k=-4:-1" = c("-4:-3","-2","-1"),
                "k=-4:-3" = c("-4:-3"),
                "k=-2:-1" = c("-2","-1"))
if (!REFERENCE %in% names(REFSETS)) stop("REFERENCE must be one of: ",
                                         paste(names(REFSETS), collapse = ", "))
REFSET <- REFSETS[[REFERENCE]]

rbin <- function(k) fifelse(k <= -11, "<=-11", fifelse(k <= -8, "-10:-8",
  fifelse(k <= -5, "-7:-5", fifelse(k <= -3, "-4:-3", fifelse(k == -2, "-2",
  fifelse(k == -1, "-1", fifelse(k == 0, "0", fifelse(k <= 3, "1:3",
  fifelse(k <= 8, "4:8", fifelse(k <= 13, "9:13",
  fifelse(k <= 18, "14:18", ">=19")))))))))))
cbin <- function(a) fifelse(a <= 2004, "2000-04", fifelse(a <= 2009, "2005-09",
  fifelse(a <= 2014, "2010-14", fifelse(a <= 2019, "2015-19", "2020-23"))))

## ------------------------------- load ---------------------------------------
message("Reading ", CROP, " field-year panel ...")
fdir <- file.path(DATA, "Field_yield_tillage", "processed_smooth_1999_2023")
fs <- list.files(fdir, pattern = paste0("^", CROP, "_yield_tillage_[A-Z]{2}_1999_2023\\.csv$"),
                 full.names = TRUE)
if (!length(fs)) stop("No panel files found in ", fdir)
if (!is.null(STATES)) {
  keep <- vapply(fs, function(f) any(vapply(STATES, function(s)
            grepl(paste0("_", s, "_"), basename(f)), logical(1))), logical(1))
  fs <- fs[keep]
  if (!length(fs)) stop("No files matched STATES = ", paste(STATES, collapse = ","))
}
message("  ", length(fs), " state file(s)")
dt <- rbindlist(lapply(fs, fread, select = c("OBJECTID","yield","size","year","till",
                                             "GEOID","till_1_count","till_1_streak")),
                use.names = TRUE)
setnames(dt, "GEOID", "FIPS")
dt[, `:=`(OBJECTID = as.character(OBJECTID), FIPS = as.integer(FIPS),
          year = as.integer(year), till = as.integer(till),
          till_1_count = as.integer(till_1_count),
          till_1_streak = as.integer(till_1_streak), size = as.numeric(size),
          log_yield = log(as.numeric(yield) / 1000))][, yield := NULL]
dt <- dt[is.finite(log_yield) & is.finite(size) & size > 0]

## rainfed only
irr <- fread(file.path(DATA, "irrigation", "field_irrigation_ever_only_irrigated.csv"))
dt <- dt[!OBJECTID %in% as.character(irr$OBJECTID)]

## drop fields already in low till at the start of the panel
omit <- unique(c(dt[year == 1999L & till_1_count == 1L, OBJECTID],
                 dt[year == 2000L & till_1_count == 2L, OBJECTID]))
dt <- dt[!OBJECTID %in% omit]

if (FIELD_FRAC < 1) {
  set.seed(SEED)
  ids <- unique(dt$OBJECTID)
  dt <- dt[OBJECTID %in% sample(ids, max(1L, floor(length(ids) * FIELD_FRAC)))]
  message("  field subsample: ", format(uniqueN(dt$OBJECTID), big.mark = ","), " fields")
}

## ------------------------- eligibility & event time -------------------------
dt[, clean := till == 0L & till_1_count == 0L]                     # never adopted
dt[, cont  := till == 1L & till_1_count == till_1_streak]          # unbroken low-till spell
st <- dt[, .(hc = any(cont), el = any(till_1_count > 0L),
             a = if (any(cont)) min(year[cont] - till_1_streak[cont] + 1L) else NA_integer_),
         by = OBJECTID]
st[, grp := fifelse(hc, "adopter", fifelse(!el, "never", "ambiguous"))]
dt <- merge(dt[clean | cont], st[, .(OBJECTID, grp, a)], by = "OBJECTID", sort = FALSE)

## weather
wx <- fread(file.path(DATA, "PRISM", "PRISM_climate_1999_2023.csv"),
            select = c("FIPS","year","GDD_4_5","ppt_4_5","GDD_6_9","ppt_6_9","EDD_6_9"))
setnames(wx, c("ppt_4_5","ppt_6_9"), c("PPT_4_5","PPT_6_9"))
wx[, `:=`(FIPS = as.integer(FIPS), year = as.integer(year))]
dt <- merge(dt, wx, by = c("FIPS","year"), all.x = TRUE, sort = FALSE)
COV <- c("GDD_4_5","PPT_4_5","GDD_6_9","PPT_6_9","EDD_6_9")
dt <- dt[Reduce(`&`, lapply(COV, function(v) is.finite(dt[[v]])))]
dt[, state := ST[as.character(FIPS %/% 1000L)]]

es <- dt[grp %in% c("adopter","never")]; rm(dt); gc()
es[, k   := fifelse(grp == "adopter", year - a, NA_integer_)]
es[, rel := fifelse(grp == "adopter", rbin(k), "NEVER")]
es[, coh := fifelse(grp == "adopter", cbin(a), "never")]
message(sprintf("Sample: %s field-years | %s fields | adopters %s | never-adopters %s (%.1f%% of fields)",
    format(nrow(es), big.mark=","), format(uniqueN(es$OBJECTID), big.mark=","),
    format(es[grp=="adopter", .N], big.mark=","), format(es[grp=="never", .N], big.mark=","),
    100 * es[grp=="never", uniqueN(OBJECTID)] / es[, uniqueN(OBJECTID)]))

## ------------------- stage 1: saturated cohort x period ---------------------
es[, cell := fifelse(grp == "adopter" & !(rel %in% REFSET), paste0(coh, "|", rel), "REF")]
cnt  <- es[cell != "REF", .(n = .N, w = sum(size)), by = .(cell, coh, rel)]
keep <- cnt[n >= MIN_CELL, cell]
if (length(keep) < nrow(cnt))
  message("  ", nrow(cnt) - length(keep), " thin cell(s) folded into the reference")
es[!(cell %in% keep), cell := "REF"]
es[, cell := relevel(factor(cell), ref = "REF")]

fe_str <- if (FE == "field_year") "OBJECTID + year" else "OBJECTID + year^state"
message("Fitting saturated model (", length(keep), " cells, FE = ", fe_str, ") ...")
m <- feols(as.formula(paste("log_yield ~ i(cell, ref='REF') +",
                            paste(COV, collapse = " + "), "|", fe_str)),
           data = es, weights = ~size, cluster = ~FIPS, lean = FALSE, mem.clean = TRUE)

## ------------- stage 2: interaction weighting (cohort area shares) ----------
b <- coef(m); V <- vcov(m)
nm  <- grep("^cell::", names(b), value = TRUE)
key <- data.table(term = nm, cell = sub("^cell::", "", nm))
key[, coh := sub("\\|.*$", "", cell)][, rel := sub("^.*\\|", "", cell)]
wts <- cnt[cell %in% keep][key, on = "cell"]

agg <- rbindlist(lapply(setdiff(RL, REFSET), function(L) {
  s <- wts[rel == L]
  if (!nrow(s)) return(NULL)
  w  <- s$w / sum(s$w)                       # cohort area shares: non-negative
  tm <- s$term
  data.table(rel = L, n_cohorts = nrow(s), estimate = sum(w * b[tm]),
             std_error = sqrt(max(as.numeric(t(w) %*% V[tm, tm, drop = FALSE] %*% w), 0)))
}))
## the omitted cells have no coefficient -- report them explicitly, pinned at zero
agg <- rbind(agg, data.table(rel = REFSET, n_cohorts = NA_integer_,
                             estimate = 0, std_error = NA_real_), fill = TRUE)
agg[, `:=`(is_reference = rel %in% REFSET,
           period   = fifelse(match(rel, RL) <= 6L, "pre", "post"),
           duration = c("0"="1","1:3"="2-4","4:8"="5-9","9:13"="10-14",
                        "14:18"="15-19",">=19"="20+")[rel],
           pct      = (exp(estimate) - 1) * 100,
           pct_lo   = (exp(estimate - 1.96 * std_error) - 1) * 100,
           pct_hi   = (exp(estimate + 1.96 * std_error) - 1) * 100)]
agg[is_reference == TRUE, `:=`(pct = 0, pct_lo = NA_real_, pct_hi = NA_real_)]
agg[, rel := factor(rel, levels = RL)]; setorder(agg, rel)
agg[, `:=`(crop = CROP, reference = REFERENCE, fe = fe_str)]

cat(sprintf("\nSun-Abraham IW event study -- %s | reference %s | FE %s | obs %s\n",
            CROP, REFERENCE, fe_str, format(m$nobs, big.mark = ",")))
cat(sprintf("  %-8s %8s %7s %11s %11s %22s\n",
            "rel k","period","cohorts","estimate","pct","95% CI (pct)"))
for (i in seq_len(nrow(agg))) {
  if (agg$is_reference[i]) {
    cat(sprintf("  %-8s %8s %7s %11s %11s %22s\n", as.character(agg$rel[i]),
                agg$period[i], "-", "reference", "0.000%", "-"))
  } else {
    cat(sprintf("  %-8s %8s %7d %+11.6f %+10.3f%% [%+8.3f, %+8.3f]\n",
        as.character(agg$rel[i]), agg$period[i], agg$n_cohorts[i], agg$estimate[i],
        agg$pct[i], agg$pct_lo[i], agg$pct_hi[i]))
  }
}
cat("\nNOTE: cells marked 'reference' are the omitted category and have no coefficient.\n")
if (all(REFSET %in% RL[1:6]) && length(REFSET) == 6L)
  cat("      This reference pins the whole pre-period, so it cannot test pre-trends.\n")
cat(sprintf("      Cohorts contributing fall with exposure (%d at the longest bin):\n",
            agg[rel == ">=19", n_cohorts]))
cat("      with a panel starting in 1999 only early adopters reach 20+ years.\n")

tag <- gsub("[^A-Za-z0-9]+", "_", paste(CROP, REFERENCE, FE, sep = "_"))
fout <- file.path(OUTDIR, paste0("sunab_", tag, ".csv"))
fwrite(agg[, .(crop, reference, fe, rel, period, duration, is_reference, n_cohorts,
               estimate, std_error, pct, pct_lo, pct_hi)], fout)
cat("\nWrote ", fout, "\n", sep = "")
