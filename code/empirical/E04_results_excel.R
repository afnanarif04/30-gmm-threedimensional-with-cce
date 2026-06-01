# =============================================================================
# E04_results_excel.R
# Empirical pipeline Step 4 — Consolidate all results into one Excel file
# Paper: [anonymous]
# Author: [anonymous]
#
# INPUT:  emp_panel.rds          (from E01)
#         cd_diagnostics.csv     (from E02)
#         estimation_results.csv (from E03)
# OUTPUT: EMPIRICAL_RESULTS.xlsx (one file, four sheets)
#
# SHEETS:
#   1.PanelSummary   — dataset description for Section 7.2
#   2.CD_Diagnostics — Table 3 top panel: CD and BKP exponents
#   3.Estimates      — Table 3 bottom panel: MG estimates with SEs
#   4.PairStats      — per-pair summary statistics (for supplement)
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr,    quietly = TRUE, warn.conflicts = FALSE)
library(openxlsx, quietly = TRUE, warn.conflicts = FALSE)

select    <- dplyr::select;  filter    <- dplyr::filter
rename    <- dplyr::rename;  mutate    <- dplyr::mutate
arrange   <- dplyr::arrange; summarise <- dplyr::summarise

# Safe reader for intermediate CSVs
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else data.frame()

# ---- Load inputs -----------------------------------------------------------
dat    <- readRDS("emp_panel.rds")
cd_res <- rd("cd_diagnostics.csv")
est    <- rd("estimation_results.csv")

panel  <- dat$panel
N <- dat$N; M <- dat$M; T_ <- dat$T_

if (nrow(cd_res) == 0) stop("cd_diagnostics.csv missing.  Run E02 first.")
if (nrow(est)    == 0) stop("estimation_results.csv missing.  Run E03 first.")

# Fix NaN stars from Python (if re-imported)
est$stars[is.na(est$stars)] <- ""

# ---- openxlsx styles -------------------------------------------------------
wb     <- createWorkbook()
hdr_st <- createStyle(textDecoration = "bold", halign = "center",
                      fgFill = "#1F3864", fontColour = "#FFFFFF",
                      border = "Bottom", borderColour = "#FFFFFF")
sub_st <- createStyle(textDecoration = "bold", halign = "left",
                      fgFill = "#D9E1F2", border = "Bottom")
note_st <- createStyle(fontSize = 9, fontColour = "#555555",
                       textDecoration = "italic")
sig_st <- createStyle(fgFill = "#D6EAD6")     # green highlight for ***

add_sheet <- function(wb, nm, data, note = NULL) {
  addWorksheet(wb, nm)
  if (nrow(data) == 0) { writeData(wb, nm, "No data."); return(invisible()) }
  writeData(wb, nm, data, startRow = 1, headerStyle = hdr_st)
  setColWidths(wb, nm, cols = seq_along(data), widths = "auto")
  if (!is.null(note)) {
    nr <- nrow(data) + 3
    writeData(wb, nm, note, startRow = nr, startCol = 1)
    addStyle(wb, nm, note_st, rows = nr, cols = 1)
  }
}

# ---- Sheet 1: Panel summary -----------------------------------------------
panel_sum <- data.frame(
  Item  = c("N (origin economies)", "M (destination economies)", "T (years)",
            "Sample period", "Balanced pairs", "Total observations",
            "Origin economies", "Destination economies",
            "Outcome variable", "Regressor 1 (X1)", "Regressor 2 (X2)",
            "Estimator", "Lag order p_T"),
  Value = c(
    as.character(N), as.character(M), as.character(T_),
    paste0(min(dat$years_available), "\u2013", max(dat$years_available)),
    as.character(length(unique(paste(panel$origin, panel$dest)))),
    as.character(nrow(panel)),
    paste(sort(unique(panel$origin)), collapse = ", "),
    paste(sort(unique(panel$dest)),   collapse = ", "),
    "ln(bilateral exports + 1), USD thousands",
    "ln(GDP, origin economy), current USD",
    "ln(GDP, destination economy), current USD",
    "Dynamic three-dimensional CCE mean group estimator",
    "0 (ratio=3.22 at T=30; threshold 2.5)"
  ),
  stringsAsFactors = FALSE
)
add_sheet(wb, "1.PanelSummary", panel_sum,
          note = paste("Sources: CEPII BACI HS92 (1995-2016) and HS17 (2017-2024), V202601.",
                       "World Bank WDI (NY.GDP.MKTP.CD). Taiwan GDP: IMF World Economic Outlook."))

# ---- Sheet 2: CD diagnostics -----------------------------------------------
cd_out <- cd_res %>%
  select(label, cd_I, pval_I, alpha_I, cd_J, pval_J, alpha_J) %>%
  rename(`Specification`    = label,
         `CD (origin)`      = cd_I,
         `p-value (origin)` = pval_I,
         `alpha_I (BKP)`    = alpha_I,
         `CD (dest)`        = cd_J,
         `p-value (dest)`   = pval_J,
         `alpha_J (BKP)`    = alpha_J)

add_sheet(wb, "2.CD_Diagnostics", cd_out,
          note = paste("CD = Pesaran (2015) cross-sectional dependence statistic.",
                       "alpha = Bailey-Kapetanios-Pesaran (2016) exponent; alpha > 0.5 indicates",
                       "strong dependence. Step 1: no augmentation. Step 2: grand average only",
                       "(Pesaran 2006). Step 3: full three-margin augmentation (proposed)."))

# ---- Sheet 3: Estimation results -------------------------------------------
# Add significance highlight for *** rows
addWorksheet(wb, "3.Estimates")
est_out <- est %>%
  select(label, estimate, se, tstat, pval, stars) %>%
  rename(`Parameter` = label, `Estimate` = estimate, `Std Error` = se,
         `t-stat` = tstat, `p-value` = pval, `Sig.` = stars)

writeData(wb, "3.Estimates", est_out, startRow = 1, headerStyle = hdr_st)
setColWidths(wb, "3.Estimates", cols = seq_along(est_out), widths = "auto")

# Highlight *** rows
for (r in seq_len(nrow(est_out))) {
  if (trimws(est_out$Sig.[r]) == "***")
    addStyle(wb, "3.Estimates", sig_st, rows = r + 1, cols = seq_along(est_out),
             gridExpand = TRUE, stack = TRUE)
}

note_row <- nrow(est_out) + 3
writeData(wb, "3.Estimates",
          paste("* p<0.10  ** p<0.05  *** p<0.01.",
                "Standard errors based on margin-based variance estimator: V = V_I + V_J.",
                "Long-run estimates via delta method.",
                sprintf("Valid pairs: %d / %d.", est$n_valid_pairs[1], N * M)),
          startRow = note_row, startCol = 1)
addStyle(wb, "3.Estimates", note_st, rows = note_row, cols = 1)

# ---- Sheet 4: Per-pair summary statistics ----------------------------------
pair_stats <- panel %>%
  group_by(origin, dest) %>%
  summarise(
    N_years       = n(),
    mean_ln_trade = round(mean(ln_trade, na.rm = TRUE), 4),
    sd_ln_trade   = round(sd(ln_trade,   na.rm = TRUE), 4),
    min_ln_trade  = round(min(ln_trade,  na.rm = TRUE), 4),
    max_ln_trade  = round(max(ln_trade,  na.rm = TRUE), 4),
    ar1_ln_trade  = round({
      y <- ln_trade[!is.na(ln_trade)]
      if (length(y) > 2 && sd(y) > 0) cor(y[-length(y)], y[-1]) else NA_real_
    }, 4),
    mean_trade_USD_bn = round(mean(trade_value, na.rm = TRUE) / 1e6, 4),
    .groups = "drop"
  ) %>%
  arrange(origin, dest)

add_sheet(wb, "4.PairStats", pair_stats,
          note = "AR(1) computed on raw ln(trade+1) series before factor adjustment.")

# ---- Save ------------------------------------------------------------------
saveWorkbook(wb, "EMPIRICAL_RESULTS.xlsx", overwrite = TRUE)
cat("Saved: EMPIRICAL_RESULTS.xlsx\n")

# Quick check
wb2 <- loadWorkbook("EMPIRICAL_RESULTS.xlsx")
cat(sprintf("Sheets: %s\n", paste(names(wb2), collapse = ", ")))
cat("\nE04_results_excel.R done.  All empirical results in EMPIRICAL_RESULTS.xlsx.\n")
cat("For Table 3 in the paper: use Sheet '3.Estimates' and '2.CD_Diagnostics'.\n")
