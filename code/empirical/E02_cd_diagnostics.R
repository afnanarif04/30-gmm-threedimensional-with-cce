# =============================================================================
# E02_cd_diagnostics.R
# Empirical pipeline Step 2 — Cross-sectional dependence diagnostics
# Paper: [anonymous]
# Author: [anonymous]
#
# INPUT:  emp_panel.rds  (from E01_data_from_excel.R)
# OUTPUT: cd_diagnostics.csv  (three-step Section 5.2 procedure)
#
# WHAT THIS SCRIPT DOES:
#   Implements the three-step hierarchical factor diagnostic of Section 5.2:
#   Step 1 — No augmentation:     CD on raw ln_trade margins
#   Step 2 — Grand average only:  Pesaran (2006) single-level augmentation
#   Step 3 — Full 3D augmentation: proposed three-margin absorber
#   Reports Pesaran (2015) CD statistic and BKP (2016) exponent per margin.
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
select    <- dplyr::select;  filter    <- dplyr::filter
rename    <- dplyr::rename;  mutate    <- dplyr::mutate
arrange   <- dplyr::arrange; summarise <- dplyr::summarise

if (!file.exists("emp_panel.rds"))
  stop("emp_panel.rds not found.  Run E01_data_from_excel.R first.")

dat <- readRDS("emp_panel.rds")
Y   <- dat$Y; X1 <- dat$X1; X2 <- dat$X2
N   <- dat$N; M  <- dat$M;  T_ <- dat$T_

cat(sprintf("Panel loaded: N=%d, M=%d, T=%d\n", N, M, T_))

# ---- Helper: Pesaran (2015) CD statistic on a K×T residual matrix ----------
pesaran_cd <- function(res_mat) {
  K <- nrow(res_mat); T <- ncol(res_mat)
  rho_sum <- 0; n_pairs <- 0L
  for (a in seq_len(K - 1)) {
    for (b in seq(a + 1, K)) {
      mask <- !is.na(res_mat[a, ]) & !is.na(res_mat[b, ])
      if (sum(mask) < 3) next
      r <- tryCatch(cor(res_mat[a, mask], res_mat[b, mask]), error = function(e) NA_real_)
      if (is.finite(r)) { rho_sum <- rho_sum + sqrt(sum(mask)) * r; n_pairs <- n_pairs + 1L }
    }
  }
  if (n_pairs == 0) return(c(cd = NA_real_, pval = NA_real_, avg_rho = NA_real_))
  cd      <- sqrt(2 / (K * (K - 1))) * rho_sum
  pval    <- 2 * pnorm(-abs(cd))
  avg_rho <- rho_sum / (n_pairs * sqrt(T))
  c(cd = round(cd, 3), pval = round(pval, 4), avg_rho = round(avg_rho, 4))
}

# ---- Helper: BKP (2016) exponent estimator ----------------------------------
bkp_exponent <- function(res_mat) {
  K  <- nrow(res_mat)
  rs <- numeric(0)
  for (a in seq_len(K - 1)) {
    for (b in seq(a + 1, K)) {
      mask <- !is.na(res_mat[a, ]) & !is.na(res_mat[b, ])
      if (sum(mask) < 3) next
      r <- tryCatch(abs(cor(res_mat[a, mask], res_mat[b, mask])),
                    error = function(e) NA_real_)
      if (is.finite(r) && r > 0) rs <- c(rs, r)
    }
  }
  if (length(rs) < 3) return(NA_real_)
  alpha <- 1 + log(mean(rs)) / log(K)
  round(max(0, min(1, alpha)), 3)
}

# ---- Margin means helper ----------------------------------------------------
make_margins <- function(Y_, X1_, X2_) {
  N_ <- dim(Y_)[1]; M_ <- dim(Y_)[2]; T_ <- dim(Y_)[3]
  list(
    ga_Y  = apply(Y_,  3, mean, na.rm = TRUE),        # [T]
    ga_X1 = apply(X1_, 3, mean, na.rm = TRUE),
    ga_X2 = apply(X2_, 3, mean, na.rm = TRUE),
    oa_Y  = apply(Y_,  c(1, 3), mean, na.rm = TRUE),  # [N x T]
    oa_X1 = apply(X1_, c(1, 3), mean, na.rm = TRUE),
    oa_X2 = apply(X2_, c(1, 3), mean, na.rm = TRUE),
    da_Y  = apply(Y_,  c(2, 3), mean, na.rm = TRUE),  # [M x T]
    da_X1 = apply(X1_, c(2, 3), mean, na.rm = TRUE),
    da_X2 = apply(X2_, c(2, 3), mean, na.rm = TRUE)
  )
}

# ---- Margin residual matrices -----------------------------------------------
origin_margins <- function(res3) apply(res3, c(1, 3), mean, na.rm = TRUE)  # [N x T]
dest_margins   <- function(res3) apply(res3, c(2, 3), mean, na.rm = TRUE)  # [M x T]

# ---- Demean functions -------------------------------------------------------
demean_grand <- function(arr) {
  mg <- make_margins(arr, X1, X2)
  out <- arr
  for (i in seq_len(N)) for (j in seq_len(M)) {
    y    <- arr[i, j, ]
    mask <- !is.na(y)
    if (sum(mask) < 5) next
    A  <- cbind(1, mg$ga_Y[mask], mg$ga_X1[mask], mg$ga_X2[mask])
    cf <- tryCatch(.lm.fit(A, y[mask])$coefficients, error = function(e) NULL)
    if (!is.null(cf) && all(is.finite(cf)))
      out[i, j, ][mask] <- y[mask] - A %*% cf
  }
  out
}

demean_full <- function(arr) {
  mg <- make_margins(arr, X1, X2)
  out <- arr
  for (i in seq_len(N)) for (j in seq_len(M)) {
    y    <- arr[i, j, ]
    mask <- !is.na(y)
    if (sum(mask) < 12) next
    A <- cbind(1,
               mg$ga_Y[mask],   mg$ga_X1[mask],   mg$ga_X2[mask],
               mg$oa_Y[i, mask], mg$oa_X1[i, mask], mg$oa_X2[i, mask],
               mg$da_Y[j, mask], mg$da_X1[j, mask], mg$da_X2[j, mask])
    cf <- tryCatch(.lm.fit(A, y[mask])$coefficients, error = function(e) NULL)
    if (!is.null(cf) && all(is.finite(cf)))
      out[i, j, ][mask] <- y[mask] - A %*% cf
  }
  out
}

# ---- Step 1: raw ln_trade ---------------------------------------------------
cat("\nStep 1 — No augmentation (raw ln_trade)...\n")
res_I0 <- origin_margins(Y)
res_J0 <- dest_margins(Y)
cd0_I  <- pesaran_cd(res_I0);  al0_I <- bkp_exponent(res_I0)
cd0_J  <- pesaran_cd(res_J0);  al0_J <- bkp_exponent(res_J0)
cat(sprintf("  Origin: CD=%.3f (p=%.4f), alpha=%.3f\n", cd0_I["cd"], cd0_I["pval"], al0_I))
cat(sprintf("  Dest:   CD=%.3f (p=%.4f), alpha=%.3f\n", cd0_J["cd"], cd0_J["pval"], al0_J))

# ---- Step 2: grand average only --------------------------------------------
cat("Step 2 — Grand average only (Pesaran 2006 baseline)...\n")
res1_Y <- demean_grand(Y)
res_I1 <- origin_margins(res1_Y);  res_J1 <- dest_margins(res1_Y)
cd1_I  <- pesaran_cd(res_I1);  al1_I <- bkp_exponent(res_I1)
cd1_J  <- pesaran_cd(res_J1);  al1_J <- bkp_exponent(res_J1)
cat(sprintf("  Origin: CD=%.3f (p=%.4f), alpha=%.3f\n", cd1_I["cd"], cd1_I["pval"], al1_I))
cat(sprintf("  Dest:   CD=%.3f (p=%.4f), alpha=%.3f\n", cd1_J["cd"], cd1_J["pval"], al1_J))

# ---- Step 3: full 3D augmentation (proposed) --------------------------------
cat("Step 3 — Full 3D augmentation (proposed)...\n")
res2_Y <- demean_full(Y)
res_I2 <- origin_margins(res2_Y);  res_J2 <- dest_margins(res2_Y)
cd2_I  <- pesaran_cd(res_I2);  al2_I <- bkp_exponent(res_I2)
cd2_J  <- pesaran_cd(res_J2);  al2_J <- bkp_exponent(res_J2)
cat(sprintf("  Origin: CD=%.3f (p=%.4f), alpha=%.3f\n", cd2_I["cd"], cd2_I["pval"], al2_I))
cat(sprintf("  Dest:   CD=%.3f (p=%.4f), alpha=%.3f\n", cd2_J["cd"], cd2_J["pval"], al2_J))

# ---- Compile and save -------------------------------------------------------
cd_diag <- data.frame(
  spec    = c("spec0", "spec1", "spec2"),
  label   = c("No augmentation", "Grand average only", "Full 3D augmentation (proposed)"),
  cd_I    = c(cd0_I["cd"],   cd1_I["cd"],   cd2_I["cd"]),
  pval_I  = c(cd0_I["pval"], cd1_I["pval"], cd2_I["pval"]),
  alpha_I = c(al0_I, al1_I, al2_I),
  cd_J    = c(cd0_J["cd"],   cd1_J["cd"],   cd2_J["cd"]),
  pval_J  = c(cd0_J["pval"], cd1_J["pval"], cd2_J["pval"]),
  alpha_J = c(al0_J, al1_J, al2_J),
  stringsAsFactors = FALSE
)

cat("\n=== DIAGNOSTICS SUMMARY ===\n")
print(cd_diag[, c("label","cd_I","alpha_I","cd_J","alpha_J")], row.names = FALSE)

write.csv(cd_diag, "cd_diagnostics.csv", row.names = FALSE)
cat("\nSaved: cd_diagnostics.csv\n")
cat("E02_cd_diagnostics.R done.  Run E03_estimation.R next.\n")
