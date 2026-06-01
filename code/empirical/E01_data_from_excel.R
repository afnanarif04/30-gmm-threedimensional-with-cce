# =============================================================================
# E01_data_from_excel.R
# Empirical pipeline Step 1 — Load and prepare data from the extracted Excel
# Paper: [anonymous]
# Author: [anonymous]
#
# INPUT:  Empirical_Data_3DCCE_Paper.xlsx  (produced by the data-assembly stage)
# OUTPUT: emp_panel.rds  — clean data object with 3D arrays for the estimator
#
# WHAT THIS SCRIPT DOES:
#   Reads Sheet "1.BilateralPanel" (trade) and Sheet "2.GDP" from the single
#   extracted Excel workbook.  No raw BACI files are needed here.
#   Builds the N x M x T arrays Y (ln trade), X1 (ln GDP origin),
#   X2 (ln GDP destination) and saves them together with the panel index tables.
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr,     quietly = TRUE, warn.conflicts = FALSE)
library(tidyr,     quietly = TRUE, warn.conflicts = FALSE)
library(openxlsx,  quietly = TRUE, warn.conflicts = FALSE)

select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
arrange   <- dplyr::arrange
summarise <- dplyr::summarise

EXCEL_FILE <- "Empirical_Data_3DCCE_Paper.xlsx"
if (!file.exists(EXCEL_FILE))
  stop(EXCEL_FILE, " not found.  Place it in the same folder as this script.")

# ---- A. Read bilateral trade panel -----------------------------------------
cat("Reading Sheet '1.BilateralPanel'...\n")
panel_raw <- read.xlsx(EXCEL_FILE, sheet = "1.BilateralPanel")
cat(sprintf("  Rows: %d,  Columns: %s\n",
            nrow(panel_raw), paste(names(panel_raw), collapse = ", ")))

# Standardise column names (the Excel header contains spaces)
names(panel_raw) <- make.names(names(panel_raw))
# Expected: Origin..ISO3., Destination..ISO3., Year, i..origin.index.,
#           j..dest.index., t..time.index., Trade.value..USD.thousands.,
#           ln.Trade.1., ln.GDP_origin., ln.GDP_dest.
panel <- panel_raw %>%
  rename(
    origin      = Origin..ISO3.,
    dest        = Destination..ISO3.,
    year        = Year,
    i           = i..origin.index.,
    j           = j..dest.index.,
    t           = t..time.index.,
    trade_value = Trade.value..USD.thousands.,
    ln_trade    = ln.Trade.1.,
    ln_gdp_o    = ln.GDP_origin.,
    ln_gdp_d    = ln.GDP_dest.
  ) %>%
  mutate(across(c(i, j, t, year), as.integer)) %>%
  arrange(i, j, t)

cat(sprintf("  Origins (N=%d): %s\n",
            length(unique(panel$origin)), paste(sort(unique(panel$origin)), collapse = ", ")))
cat(sprintf("  Dests   (M=%d): %s\n",
            length(unique(panel$dest)),   paste(sort(unique(panel$dest)), collapse = ", ")))
cat(sprintf("  Years   (T=%d): %d-%d\n",
            length(unique(panel$year)),   min(panel$year), max(panel$year)))
cat(sprintf("  Obs: %d,  Balanced pairs: %d\n",
            nrow(panel), length(unique(paste(panel$origin, panel$dest)))))

# ---- B. Read GDP series (for reference / completeness check) ---------------
cat("\nReading Sheet '2.GDP'...\n")
gdp_raw <- read.xlsx(EXCEL_FILE, sheet = "2.GDP")
names(gdp_raw) <- make.names(names(gdp_raw))
gdp <- gdp_raw %>%
  rename(iso3c = Country..ISO3., year = Year, gdp = GDP..current.USD., ln_gdp = ln.GDP.) %>%
  mutate(year = as.integer(year)) %>%
  filter(!is.na(ln_gdp))
cat(sprintf("  GDP obs: %d,  Countries: %d,  Years: %d-%d\n",
            nrow(gdp), length(unique(gdp$iso3c)), min(gdp$year), max(gdp$year)))

# ---- C. Index tables --------------------------------------------------------
N      <- max(panel$i)
M      <- max(panel$j)
T_     <- max(panel$t)
years  <- sort(unique(panel$year))

origins_idx <- panel %>%
  distinct(origin, i) %>% arrange(i)
dests_idx   <- panel %>%
  distinct(dest,   j) %>% arrange(j)
years_idx   <- panel %>%
  distinct(year,   t) %>% arrange(t)

cat(sprintf("\nPanel dimensions: N=%d, M=%d, T=%d\n", N, M, T_))

# ---- D. Build 3D arrays: Y[N, M, T], X1[N, M, T], X2[N, M, T] -------------
cat("Building 3D arrays...\n")
Y  <- array(NA_real_, dim = c(N, M, T_))
X1 <- array(NA_real_, dim = c(N, M, T_))
X2 <- array(NA_real_, dim = c(N, M, T_))

for (row_k in seq_len(nrow(panel))) {
  ii <- panel$i[row_k]
  jj <- panel$j[row_k]
  tt <- panel$t[row_k]
  Y [ii, jj, tt] <- panel$ln_trade[row_k]
  X1[ii, jj, tt] <- panel$ln_gdp_o[row_k]
  X2[ii, jj, tt] <- panel$ln_gdp_d[row_k]
}

nan_Y  <- sum(is.na(Y))
nan_X1 <- sum(is.na(X1))
cat(sprintf("  Y: %s,  NaN=%d  (%.1f%%)\n",
            paste(dim(Y), collapse = " x "), nan_Y, 100 * nan_Y / prod(dim(Y))))
cat(sprintf("  X1: %s, NaN=%d\n", paste(dim(X1), collapse = " x "), nan_X1))

# ---- E. Persistence check (confirms dynamic specification is warranted) -----
cat("\nPersistence check (pair-level AR(1) of ln_trade):\n")
ar1_vec <- numeric(0)
for (ii in seq_len(N)) {
  for (jj in seq_len(M)) {
    y_ij <- Y[ii, jj, ]
    valid <- !is.na(y_ij)
    if (sum(valid) < 3 || sd(y_ij[valid]) < 1e-8) next
    yy <- y_ij[valid]
    r  <- tryCatch(cor(yy[-length(yy)], yy[-1]), error = function(e) NA_real_)
    if (is.finite(r)) ar1_vec <- c(ar1_vec, r)
  }
}
cat(sprintf("  Median AR(1): %.3f,  Share > 0.80: %.1f%%\n",
            median(ar1_vec), 100 * mean(ar1_vec > 0.80)))
cat("  (Raw AR(1) includes common-shock persistence; factor-adjusted rho_MG will be lower)\n")

# ---- F. Save ----------------------------------------------------------------
saveRDS(
  list(Y = Y, X1 = X1, X2 = X2,
       panel       = panel,
       gdp         = gdp,
       origins_idx = origins_idx,
       dests_idx   = dests_idx,
       years_idx   = years_idx,
       N = N, M = M, T_ = T_,
       years_available = years),
  "emp_panel.rds"
)
cat("\nSaved: emp_panel.rds\n")
cat("E01_data_from_excel.R done.  Run E02_cd_diagnostics.R next.\n")
