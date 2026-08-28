###############################################################################
##  Cumulative low-till effect on log yield -- CORN
##  Companion to 1_1 (continuous streak). Two-way fixed effects, field and year.
##
##  The published models use till_1_streak on a restricted sample: eligibility is
##  clean_control | continuous_low_till, and the second condition requires
##  till_1_count == till_1_streak, i.e. the field has NEVER reverted. That drops
##  roughly half the panel. These models use till_1_count (cumulative low-till
##  years, ever) and keep the reverted field-years.
##
##  Seven specifications, mirroring the paper's linear/bin pair:
##    (1) published : streak,      restricted   -- benchmark
##    (2) identity   : cumulative, restricted   -- MUST equal (1); see note below
##    (3) cumulative : cumulative, full sample
##    (4) stock/flow : cumulative + currently-low-till indicator, full sample
##    (5) published : streak bins, restricted   -- benchmark
##    (6) cumulative bins, full sample
##    (7) cumulative bins + currently-low-till, full sample
##
##  NOTE on (2). On the restricted sample the eligibility rule forces
##  till_1_count == till_1_streak row by row, so cumulative and streak are the SAME
##  column and (2) reproduces (1) to machine precision. It is included precisely to
##  prove that everything which changes between (1) and (3) comes from the SAMPLE,
##  not from the treatment definition.
###############################################################################

## =============================== PORTABLE PATHS ==============================
## Set ROOT_OVERRIDE to the folder that contains "data/" and leave everything else
## alone. If you leave it as "", the script locates the project automatically.
##
## Required layout under that folder:
##   data/Field_yield_tillage/processed_smooth_1999_2023/corn_yield_tillage_XX_1999_2023.csv
##   data/PRISM/PRISM_climate_1999_2023.csv
##   data/irrigation/field_irrigation_ever_only_irrigated.csv
## ============================================================================
ROOT_OVERRIDE <- ""

CROP       <- "corn"
STATES     <- NULL      # NULL = all nine states; else e.g. c("IA","IL")
FIELD_FRAC <- 1         # 1 = all fields; else a fraction, e.g. 0.25
SEED       <- 42
## MEMORY: the full-sample models run on ~44M corn field-years and peak near
## 30-50 GB. On a laptop set STATES <- c("IA") and/or FIELD_FRAC <- 0.25.

pkgs <- c("data.table", "fixest")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing package(s): ", paste(miss, collapse = ", "),
                       "\nInstall with: install.packages(c(",
                       paste(sprintf('"%s"', miss), collapse = ", "), "))")
library(data.table); library(fixest)
setDTthreads(0L)

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
  looks_like_root <- function(p) {
    dir.exists(file.path(p, "data", "Field_yield_tillage",
                         "processed_smooth_1999_2023")) &&
      file.exists(file.path(p, "data", "PRISM", "PRISM_climate_1999_2023.csv"))
  }
  candidates <- character()
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
  env <- Sys.getenv("PROJECT_ROOT", unset = "")
  if (nzchar(env)) candidates <- c(candidates, ascend_paths(env))
  frame_files <- character()
  for (fr in sys.frames()) {
    of <- tryCatch(fr$ofile, error = function(e) NULL)
    if (!is.null(of)) frame_files <- c(frame_files, as.character(of))
  }
  if (!length(frame_files) && requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (is.character(p) && length(p) == 1L && nzchar(p)) frame_files <- p
  }
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(fa)) frame_files <- c(frame_files, sub("^--file=", "", fa[1]))
  for (f in unique(frame_files)) candidates <- c(candidates, ascend_paths(dirname(f)))
  candidates <- unique(c(candidates, ascend_paths(getwd())))
  for (cand in candidates) if (looks_like_root(cand)) return(normalizePath(cand))
  stop("Could not locate the project root.\n",
       "Set ROOT_OVERRIDE at the top of this script to the folder containing data/.\n",
       "Searched ", length(candidates), " candidate folder(s), starting from:\n  ",
       paste(utils::head(candidates, 4), collapse = "\n  "))
}

root_dir   <- find_root_dir()
data_dir   <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "output_linear_1999_2023",
                        paste0("4_", CROP, "_cumulative_field_year"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Project root : ", root_dir)
message("Output folder: ", output_dir)

COV  <- c("GDD_4_5", "PPT_4_5", "GDD_6_9", "PPT_6_9", "EDD_6_9")
COVf <- paste(COV, collapse = " + ")
BL   <- c("0", "1_4", "5_9", "10_14", "15_19", "20+")
binit <- function(x) factor(fifelse(x == 0, "0", fifelse(x <= 4, "1_4",
  fifelse(x <= 9, "5_9", fifelse(x <= 14, "10_14",
  fifelse(x <= 19, "15_19", "20+"))))), levels = BL)

## ------------------------------- load ---------------------------------------
fdir <- file.path(data_dir, "Field_yield_tillage", "processed_smooth_1999_2023")
fs <- list.files(fdir, pattern = paste0("^", CROP, "_yield_tillage_[A-Z]{2}_1999_2023\\.csv$"),
                 full.names = TRUE)
if (!length(fs)) stop("No panel files found in ", fdir)
if (!is.null(STATES)) {
  keep <- vapply(fs, function(f) any(vapply(STATES, function(s)
            grepl(paste0("_", s, "_"), basename(f)), logical(1))), logical(1))
  fs <- fs[keep]
  if (!length(fs)) stop("No files matched STATES = ", paste(STATES, collapse = ","))
}
message("Reading ", length(fs), " state file(s) for ", CROP, " ...")
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
irr <- fread(file.path(data_dir, "irrigation", "field_irrigation_ever_only_irrigated.csv"))
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

## weather
wx <- fread(file.path(data_dir, "PRISM", "PRISM_climate_1999_2023.csv"),
            select = c("FIPS","year","GDD_4_5","ppt_4_5","GDD_6_9","ppt_6_9","EDD_6_9"))
setnames(wx, c("ppt_4_5","ppt_6_9"), c("PPT_4_5","PPT_6_9"))
wx[, `:=`(FIPS = as.integer(FIPS), year = as.integer(year))]
dt <- merge(dt, wx, by = c("FIPS","year"), all.x = TRUE, sort = FALSE)
dt <- dt[Reduce(`&`, lapply(COV, function(v) is.finite(dt[[v]])))]

## ---------------------------- treatment variables ---------------------------
dt[, clean   := till == 0L & till_1_count == 0L]                  # never adopted
dt[, cont    := till == 1L & till_1_count == till_1_streak]        # unbroken spell
dt[, cum     := as.numeric(till_1_count)]                          # cumulative years
dt[, streak  := fifelse(cont, as.numeric(till_1_streak), 0)]       # published treatment
dt[, current := as.numeric(till == 1L)]                            # currently low till
dt[, cum_bin    := binit(till_1_count)]
dt[, streak_bin := binit(fifelse(cont, till_1_streak, 0L))]
sub <- dt[clean | cont]                                            # published sample

## sample composition, for the SI
comp <- dt[, .(observations = .N, fields = uniqueN(OBJECTID),
               acres_M = round(sum(size) / 1e6, 1),
               mean_cumulative = round(mean(cum), 2)),
           by = .(status = fifelse(till_1_count == 0L, "never adopted",
                            fifelse(till == 1L, "currently low till",
                                    "reverted (currently high till)")))]
setorder(comp, -observations)
cat("\nSample composition (full sample)\n"); print(comp)
cat(sprintf("published (restricted) sample: %s of %s field-years (%.1f%%)\n",
    format(nrow(sub), big.mark = ","), format(nrow(dt), big.mark = ","),
    100 * nrow(sub) / nrow(dt)))
fwrite(comp, file.path(output_dir, paste0(CROP, "_cumulative_sample_composition.csv")))

## ------------------------------- estimation ---------------------------------
F1 <- function(rhs, d) feols(as.formula(paste("log_yield ~", rhs, "+", COVf,
        "| OBJECTID + year")), data = d, weights = ~size, cluster = ~FIPS,
        mem.clean = TRUE)

grab <- function(m, label, sample_lab, focal) {
  ct <- as.data.table(coeftable(m), keep.rownames = "term")
  setnames(ct, c("term","estimate","std_error","t_value","p_value"))
  ct[, `:=`(model = label, sample = sample_lab, crop = CROP, nobs = m$nobs,
            adj_r2 = fitstat(m, "ar2")$ar2, within_r2 = fitstat(m, "wr2")$wr2,
            pct_effect = (exp(estimate) - 1) * 100,
            pct_lo = (exp(estimate - 1.96 * std_error) - 1) * 100,
            pct_hi = (exp(estimate + 1.96 * std_error) - 1) * 100)]
  cat(sprintf("\n%s  [%s]  obs %s | adj R2 %.6f | within R2 %.6f\n", label, sample_lab,
              format(m$nobs, big.mark = ","), ct$adj_r2[1], ct$within_r2[1]))
  for (v in focal) {
    r <- ct[term == v]
    if (nrow(r)) cat(sprintf("   %-18s %+11.6f (%.6f)  %+8.4f%%  [%+7.4f, %+7.4f]  p=%.2e\n",
        v, r$estimate, r$std_error, r$pct_effect, r$pct_lo, r$pct_hi, r$p_value))
  }
  ct
}

res <- list()
res[[1]] <- grab(F1("streak", sub), "(1) published: continuous streak", "restricted", "streak")
res[[2]] <- grab(F1("cum", sub),    "(2) identity check: cumulative",   "restricted", "cum")
res[[3]] <- grab(F1("cum", dt),     "(3) cumulative",                   "full", "cum")
res[[4]] <- grab(F1("cum + current", dt), "(4) cumulative + current practice",
                 "full", c("cum","current"))

## (2) must reproduce (1) exactly -- verify rather than assert
b1 <- res[[1]][term == "streak", estimate]; b2 <- res[[2]][term == "cum", estimate]
cat(sprintf("\nIdentity check (1) vs (2): %.12f vs %.12f -> %s\n", b1, b2,
    if (isTRUE(all.equal(b1, b2, tolerance = 1e-12)))
      "IDENTICAL, so all later differences come from the SAMPLE" else "DIFFER (investigate)"))

bins <- function(m, label, sample_lab) {
  ct <- grab(m, label, sample_lab, character(0))
  for (b in BL[-1]) {
    ## exact suffix match: "+" in "20+" is a regex quantifier, so grep() would
    ## silently skip that bin
    v <- ct$term[endsWith(ct$term, paste0("::", b))]
    if (length(v)) { r <- ct[term == v[1]]
      cat(sprintf("   %-12s %+8.4f%%  [%+7.4f, %+7.4f]  p=%.2e\n", b,
          r$pct_effect, r$pct_lo, r$pct_hi, r$p_value)) }
  }
  v <- ct[term == "current"]
  if (nrow(v)) cat(sprintf("   %-12s %+8.4f%%  [%+7.4f, %+7.4f]  p=%.2e\n", "current",
      v$pct_effect, v$pct_lo, v$pct_hi, v$p_value))
  ct
}
res[[5]] <- bins(F1("i(streak_bin, ref='0')", sub), "(5) published: streak bins", "restricted")
res[[6]] <- bins(F1("i(cum_bin, ref='0')", dt),     "(6) cumulative bins", "full")
res[[7]] <- bins(F1("i(cum_bin, ref='0') + current", dt),
                 "(7) cumulative bins + current practice", "full")

out <- rbindlist(res, use.names = TRUE, fill = TRUE)
setcolorder(out, c("crop","model","sample","term","estimate","std_error","t_value",
                   "p_value","pct_effect","pct_lo","pct_hi","nobs","adj_r2","within_r2"))
fout <- file.path(output_dir, paste0(CROP, "_cumulative_terms.csv"))
fwrite(out, fout)
cat("\nWrote ", fout, "\n", sep = "")
cat("Wrote ", file.path(output_dir, paste0(CROP, "_cumulative_sample_composition.csv")),
    "\n", sep = "")
