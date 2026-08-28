###############################################################################
##  Leave-one-state-out validation of the low-till duration models
##
##  Refits the published linear model AND the duration-bin model with field and
##  year fixed effects, dropping all fields from one state at a time. This is the
##  robustness check behind Fig. S5: it asks whether the duration effect is driven
##  by any single region.
##
##  For each crop the script fits 10 samples (full + 9 leave-one-out) x 2 models
##  = 20 regressions, so 40 in total.
##
##  It also reports the check that matters for the claim in the text: whether each
##  leave-one-out confidence interval contains the full-sample estimate. If every
##  interval does, no single state overturns the result.
##
##  ---------------------------------------------------------------------------
##  MEMORY / RUNTIME
##  ---------------------------------------------------------------------------
##  The full nine-state run needs roughly 30-50 GB and takes on the order of an
##  hour for both crops. On a smaller machine set FIELD_FRAC below 1. Restricting
##  STATES also works but changes the meaning of the exercise -- you are then
##  leaving one state out of a smaller universe.
###############################################################################

## =============================== PORTABLE PATHS ==============================
## Set ROOT_OVERRIDE to the folder that contains "data/" and leave everything else
## alone. If you leave it as "", the script locates the project automatically.
##
##
## Required layout under that folder:
##   data/Field_yield_tillage/processed_smooth_1999_2023/<crop>_yield_tillage_XX_1999_2023.csv
##   data/PRISM/PRISM_climate_1999_2023.csv
##   data/irrigation/field_irrigation_ever_only_irrigated.csv
## ============================================================================
ROOT_OVERRIDE <- ""

CROPS      <- c("corn", "soy")
STATES     <- NULL      # NULL = all nine. A subset still works but see the note above.
FIELD_FRAC <- 1         # 1 = all fields; else a fraction, e.g. 0.25
SEED       <- 42
MAKE_FIGURE <- TRUE     # needs ggplot2 >= 3.4 (linewidth=); skipped if not installed

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
                        "8_leave_one_state_out_field_year")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
message("Project root : ", root_dir)
message("Output folder: ", output_dir)

COV  <- c("GDD_4_5", "PPT_4_5", "GDD_6_9", "PPT_6_9", "EDD_6_9")
COVf <- paste(COV, collapse = " + ")
BL   <- c("0", "1_4", "5_9", "10_14", "15_19", "20+")
BLAB <- c("0", "1-4", "5-9", "10-14", "15-19", "20+")
binit <- function(x) factor(fifelse(x == 0, "0", fifelse(x <= 4, "1_4",
  fifelse(x <= 9, "5_9", fifelse(x <= 14, "10_14",
  fifelse(x <= 19, "15_19", "20+"))))), levels = BL)
ST <- c("17"="IL","18"="IN","19"="IA","26"="MI","27"="MN",
        "29"="MO","39"="OH","46"="SD","55"="WI")

wx <- fread(file.path(data_dir, "PRISM", "PRISM_climate_1999_2023.csv"),
            select = c("FIPS","year","GDD_4_5","ppt_4_5","GDD_6_9","ppt_6_9","EDD_6_9"))
setnames(wx, c("ppt_4_5","ppt_6_9"), c("PPT_4_5","PPT_6_9"))
wx[, `:=`(FIPS = as.integer(FIPS), year = as.integer(year))]
irr_ids <- as.character(fread(file.path(data_dir, "irrigation",
              "field_irrigation_ever_only_irrigated.csv"))$OBJECTID)

build <- function(crop) {
  fdir <- file.path(data_dir, "Field_yield_tillage", "processed_smooth_1999_2023")
  fs <- list.files(fdir, pattern = paste0("^", crop,
        "_yield_tillage_[A-Z]{2}_1999_2023\\.csv$"), full.names = TRUE)
  if (!length(fs)) stop("No panel files found in ", fdir)
  if (!is.null(STATES)) {
    keep <- vapply(fs, function(f) any(vapply(STATES, function(s)
              grepl(paste0("_", s, "_"), basename(f)), logical(1))), logical(1))
    fs <- fs[keep]
    if (!length(fs)) stop("No files matched STATES = ", paste(STATES, collapse = ","))
  }
  message("Reading ", length(fs), " state file(s) for ", crop, " ...")
  dt <- rbindlist(lapply(fs, fread, select = c("OBJECTID","yield","size","year","till",
                       "GEOID","till_1_count","till_1_streak")), use.names = TRUE)
  setnames(dt, "GEOID", "FIPS")
  dt[, `:=`(OBJECTID = as.character(OBJECTID), FIPS = as.integer(FIPS),
            year = as.integer(year), till = as.integer(till),
            till_1_count = as.integer(till_1_count),
            till_1_streak = as.integer(till_1_streak), size = as.numeric(size),
            log_yield = log(as.numeric(yield) / 1000))][, yield := NULL]
  dt <- dt[is.finite(log_yield) & is.finite(size) & size > 0]
  dt <- dt[!OBJECTID %in% irr_ids]                                   # rainfed only
  omit <- unique(c(dt[year == 1999L & till_1_count == 1L, OBJECTID], # early adopters
                   dt[year == 2000L & till_1_count == 2L, OBJECTID]))
  dt <- dt[!OBJECTID %in% omit]
  if (FIELD_FRAC < 1) {
    set.seed(SEED)
    ids <- unique(dt$OBJECTID)
    dt <- dt[OBJECTID %in% sample(ids, max(1L, floor(length(ids) * FIELD_FRAC)))]
    message("  field subsample: ", format(uniqueN(dt$OBJECTID), big.mark = ","), " fields")
  }
  dt[, clean := till == 0L & till_1_count == 0L]
  dt[, cont  := till == 1L & till_1_count == till_1_streak]
  dt <- dt[clean | cont]
  dt[, dur  := as.numeric(fifelse(clean, 0L, till_1_streak))]
  dt[, dbin := binit(fifelse(clean, 0L, till_1_streak))]
  dt[, state := ST[as.character(FIPS %/% 1000L)]]
  dt <- merge(dt, wx, by = c("FIPS","year"), all.x = TRUE, sort = FALSE)
  dt[Reduce(`&`, lapply(COV, function(v) is.finite(dt[[v]])))]
}

pctf <- function(b) (exp(b) - 1) * 100
lin_res <- list(); bin_res <- list()

for (crop in CROPS) {
  cat("\n", strrep("=", 96), "\n  CROP = ", crop, "\n", strrep("=", 96), "\n", sep = "")
  dt <- build(crop)
  states <- sort(unique(dt$state))
  cat(sprintf("full sample: %s field-years | %s fields | %d states\n\n",
      format(nrow(dt), big.mark = ","), format(uniqueN(dt$OBJECTID), big.mark = ","),
      length(states)))
  cat(sprintf("  %-8s %12s %13s | %s\n", "dropped", "obs", "linear %/yr",
              paste(sprintf("%9s", BLAB[-1]), collapse = "")))

  for (s in c("none", states)) {
    d <- if (s == "none") dt else dt[state != s]

    m1 <- feols(as.formula(paste("log_yield ~ dur +", COVf, "| OBJECTID + year")),
                data = d, weights = ~size, cluster = ~FIPS, lean = TRUE, mem.clean = TRUE)
    c1 <- coeftable(m1)
    lin_res[[paste(crop, s)]] <- data.table(
      crop = crop, dropped = s, model = "linear", nobs = m1$nobs,
      b = c1["dur", 1], se = c1["dur", 2], p = c1["dur", 4],
      pct = pctf(c1["dur", 1]),
      pct_lo = pctf(c1["dur", 1] - 1.96 * c1["dur", 2]),
      pct_hi = pctf(c1["dur", 1] + 1.96 * c1["dur", 2]))

    m2 <- feols(as.formula(paste("log_yield ~ i(dbin, ref='0') +", COVf,
                "| OBJECTID + year")), data = d, weights = ~size, cluster = ~FIPS,
                lean = TRUE, mem.clean = TRUE)
    c2 <- coeftable(m2); line <- ""
    for (b in BL[-1]) {
      ## Exact name match, on purpose. Do NOT rewrite this as grep()/regex: the
      ## "+" in "20+" is a quantifier, so grep("dbin::20+", ...) silently misses
      ## the 20+ bin and it disappears from the output without any warning.
      v <- paste0("dbin::", b)
      if (v %in% rownames(c2)) {
        bb <- c2[v, 1]; ss <- c2[v, 2]
        bin_res[[paste(crop, s, b)]] <- data.table(
          crop = crop, dropped = s, model = "bins", bin = b, nobs = m2$nobs,
          b = bb, se = ss, p = c2[v, 4], pct = pctf(bb),
          pct_lo = pctf(bb - 1.96 * ss), pct_hi = pctf(bb + 1.96 * ss))
        line <- paste0(line, sprintf("%8.3f%s", pctf(bb),
                       ifelse(c2[v, 4] < 0.05, "*", " ")))
      } else line <- paste0(line, sprintf("%9s", "n/e"))
    }
    cat(sprintf("  %-8s %12s %+12.4f  | %s\n", s, format(m1$nobs, big.mark = ","),
                pctf(c1["dur", 1]), line))
    rm(m1, m2, c1, c2); if (s != "none") rm(d); gc()
  }
  rm(dt); gc()
}

L <- rbindlist(lin_res, use.names = TRUE, fill = TRUE)
B <- rbindlist(bin_res, use.names = TRUE, fill = TRUE)
## Filenames and column names (b, se, p, pct, pct_lo, pct_hi) match the CSVs the
## revision figures were built from, so FigA_loo_linear.ipynb reads these as-is.
fwrite(L, file.path(output_dir, "8_loo_linear.csv"))
fwrite(B, file.path(output_dir, "8_loo_bins.csv"))
fwrite(rbindlist(list(L, B), use.names = TRUE, fill = TRUE),
       file.path(output_dir, "8_loo_all_estimates.csv"))

## ---- the robustness check the manuscript text relies on ---------------------
cat("\n", strrep("=", 96), "\n  Does any leave-one-out CI EXCLUDE the full-sample estimate?\n",
    strrep("=", 96), "\n", sep = "")
full <- L[dropped == "none", .(crop, ref = pct)]
chk <- merge(L[dropped != "none"], full, by = "crop")
chk[, excludes := ref < pct_lo | ref > pct_hi]
cat(sprintf("  linear model: %d of %d leave-one-out CIs exclude it\n",
            chk[excludes == TRUE, .N], nrow(chk)))
for (cr in CROPS) {
  r <- L[crop == cr & dropped != "none"]
  cat(sprintf("  %-5s full sample %+.4f %%/yr | range across drops %+.4f to %+.4f\n",
      cr, L[crop == cr & dropped == "none", pct], min(r$pct), max(r$pct)))
}
if (chk[excludes == TRUE, .N] > 0) {
  cat("\n  states whose exclusion moves the estimate outside its own CI:\n")
  print(chk[excludes == TRUE, .(crop, dropped, pct, pct_lo, pct_hi, ref)])
} else {
  cat("  -> no single state overturns the result for either crop.\n")
}
b20 <- B[bin == "20+"]
for (cr in CROPS) {
  r <- b20[crop == cr & dropped != "none"]
  cat(sprintf("  %-5s 20+ bin: full sample %+.3f%% | range across drops %+.3f to %+.3f\n",
      cr, b20[crop == cr & dropped == "none", pct], min(r$pct), max(r$pct)))
}
cat("  NOTE the bins are markedly less stable than the linear slope; report both ranges.\n")

## ---- optional figure -------------------------------------------------------
if (MAKE_FIGURE && requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
  INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTE <- "#8a8a86"; GRID <- "#e6e6e2"; GRAY <- "#b9b9b4"
  CROPLAB <- c(corn = "Corn", soy = "Soybean")
  LA <- L[dropped != "none"]
  LA[, cropf := factor(CROPLAB[crop], levels = unname(CROPLAB[CROPS]))]
  ## order rows by the size of the dropped state, so both panels share one
  ## meaningful order rather than being sorted by one crop's estimate
  fulln <- L[dropped == "none", .(crop, fulln = nobs)]
  LA <- merge(LA, fulln, by = "crop")[, dropped_obs := fulln - nobs]
  lev <- rev(LA[, .(sz = max(dropped_obs)), by = dropped][order(-sz), dropped])
  LA[, lab := factor(dropped, levels = lev)]
  fullp <- L[dropped == "none"][, cropf := factor(CROPLAB[crop], levels = unname(CROPLAB[CROPS]))]
  p <- ggplot(LA, aes(pct, lab)) +
    geom_vline(data = fullp, aes(xintercept = pct), colour = ORANGE,
               linewidth = 0.7, linetype = "22") +
    ## plain segments rather than geom_errorbarh(): the caps are set to zero width
    ## anyway, and geom_errorbarh() is deprecated as of ggplot2 4.0
    geom_segment(aes(x = pct_lo, xend = pct_hi, y = lab, yend = lab),
                 colour = GRAY, linewidth = 0.7) +
    geom_point(colour = BLUE, size = 2.2) +
    facet_wrap(~cropf, nrow = 1, scales = "free_x") +
    labs(title = "Leave-one-state-out: linear duration effect (field + year FE)",
         subtitle = "Each row drops one state, ordered by state size. Dashed orange = full-sample estimate. Bars are 95% CIs.",
         x = "Duration effect (%/yr)", y = "state dropped") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
          axis.line.x = element_line(colour = MUTE, linewidth = 0.3),
          axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
          strip.text = element_text(colour = INK, face = "bold", size = 9.5),
          strip.background = element_blank(),
          plot.title = element_text(colour = INK, face = "bold", size = 12),
          plot.subtitle = element_text(colour = INK2, size = 8.6))
  ggsave(file.path(output_dir, "FigA_loo_linear.png"), p, width = 9.6, height = 4.4, dpi = 320)
  ggsave(file.path(output_dir, "FigA_loo_linear.pdf"), p, width = 9.6, height = 4.4)

  ## Figure B: the bin profile. The nine leave-one-out curves are ONE ensemble
  ## drawn in a single neutral gray, not nine cycled hues -- the reader is meant
  ## to see the width of the envelope, not to identify individual states.
  BB <- B[, .(cropf = factor(CROPLAB[crop], levels = unname(CROPLAB[CROPS])),
              dropped, x = match(bin, BL) - 1L, pct, bin)]
  anch <- unique(BB[, .(cropf, dropped)])[, `:=`(x = 0L, pct = 0, bin = "0")]
  BB <- rbind(BB, anch, use.names = TRUE)
  setorder(BB, cropf, dropped, x)
  fullB <- BB[dropped == "none"]; looB <- BB[dropped != "none"]
  ## label only the two most deviant drops at the 20+ bin; labelling all nine
  ## would collide
  dev <- merge(looB[bin == "20+", .(cropf, dropped, pct)],
               fullB[bin == "20+", .(cropf, ref = pct)], by = "cropf")
  dev[, ad := abs(pct - ref)]; setorder(dev, cropf, -ad)
  top <- dev[, head(.SD, 2), by = cropf]
  pB <- ggplot() +
    geom_hline(yintercept = 0, colour = MUTE, linewidth = 0.35) +
    geom_line(data = looB, aes(x, pct, group = dropped), colour = GRAY, linewidth = 0.45) +
    geom_line(data = fullB, aes(x, pct), colour = BLUE, linewidth = 1.2) +
    geom_point(data = fullB, aes(x, pct), colour = BLUE, size = 1.9) +
    geom_text(data = top, aes(x = 5, y = pct, label = dropped), colour = ORANGE,
              hjust = -0.25, size = 2.9, fontface = "bold") +
    facet_wrap(~cropf, nrow = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(BL) - 1L, labels = BLAB,
                       expand = expansion(mult = c(0.04, 0.14))) +
    labs(title = "Leave-one-state-out: duration-bin profile (field + year FE)",
         subtitle = paste("Blue = full sample. Gray = one curve per dropped state.",
                          "Orange labels mark the two most deviant drops at the 20+ bin."),
         x = "Consecutive years of low till", y = "Yield effect (%)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
          axis.line.x = element_line(colour = MUTE, linewidth = 0.3),
          axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2),
          strip.text = element_text(colour = INK, face = "bold", size = 9.5),
          strip.background = element_blank(),
          plot.title = element_text(colour = INK, face = "bold", size = 12),
          plot.subtitle = element_text(colour = INK2, size = 8.6))
  ggsave(file.path(output_dir, "FigB_loo_bins.png"), pB, width = 9.6, height = 4.4, dpi = 320)
  ggsave(file.path(output_dir, "FigB_loo_bins.pdf"), pB, width = 9.6, height = 4.4)
  cat("\nWrote FigA_loo_linear and FigB_loo_bins (.png / .pdf)\n")
} else if (MAKE_FIGURE) {
  cat("\nggplot2 not installed; CSVs written, figures skipped.\n")
}

cat("\nOutputs in ", output_dir, "\n", sep = "")
print(list.files(output_dir))
cat("DONE\n")
