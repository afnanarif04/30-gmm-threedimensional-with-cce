# =============================================================================
# E03_estimation.R
# Empirical pipeline Step 3 — Dynamic 3D-CCE Mean Group Estimation
# Paper: [anonymous]
# Author: [anonymous]
#
# INPUT:  emp_panel.rds         (from E01_data_from_excel.R)
#         cd_diagnostics.csv    (from E02_cd_diagnostics.R)
# OUTPUT: estimation_results.csv
#
# WHAT THIS SCRIPT DOES:
#   1. Runs the dynamic 3D-CCE mean group estimator on the bilateral panel
#   2. Computes persistence (rho_MG), short-run (beta1, beta2), and long-run
#      (theta1, theta2) mean group estimates
#   3. Computes margin-based variance estimates + delta-method long-run SEs
#   4. Saves estimation_results.csv ready for E04_results_excel.R
#
# ECONOMETRIC NOTES:
#   - p_T = get_pT(T=30, k_reg=4, N_aug=9) = 0 (ratio=3.22, above threshold)
#   - Bivariate: k_reg=4 (intercept + rho + beta1 + beta2)
#                N_aug_per_lag=9 (Y, X1, X2 each over grand/origin/dest margins)
#   - Variance: Sigma_MG = V_I + V_J  (margin-based, Theorem 2)
#   - Long-run via delta method: theta = beta/(1-rho)
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
Y   <- dat$Y;  X1 <- dat$X1;  X2 <- dat$X2
N   <- dat$N;  M  <- dat$M;   T_ <- dat$T_

cat(sprintf("Panel loaded: N=%d, M=%d, T=%d\n", N, M, T_))

# ---- A. Lag-truncation order (bivariate) -----------------------------------
# k_reg = 4 (intercept + rho + beta1 + beta2)
# N_aug = 9 (3 variables × 3 margins per lag)
# min_ratio = 2.5 ensures pair regression is identified
get_pT <- function(T, N_aug = 9, k_reg = 4, min_ratio = 2.5) {
  p_cube <- floor(T^(1 / 3))
  p_safe <- max(0L, as.integer(
    (T - 1 - min_ratio * (k_reg + N_aug)) / (1 + N_aug * min_ratio)
  ))
  min(p_cube, p_safe)
}
p_T      <- get_pT(T_)
# t_start: 0-indexed first position of the lag (R array is 1-indexed)
# Response rows in R:  (t_start + 1) : T_
# Lag rows in R:        t_start       : (T_ - 1)
t_start  <- p_T + 1L
T_eff    <- T_ - t_start
n_params <- 4L + 9L * (p_T + 1L)
cat(sprintf("p_T=%d  T_eff=%d  n_params=%d  ratio=%.2f\n",
            p_T, T_eff, n_params, T_eff / n_params))

# ---- B. Pre-compute all cross-sectional averages ---------------------------
ga_Y  <- apply(Y,  3, mean, na.rm = TRUE)     # [T]  grand average
ga_X1 <- apply(X1, 3, mean, na.rm = TRUE)
ga_X2 <- apply(X2, 3, mean, na.rm = TRUE)

oa_Y  <- apply(Y,  c(1, 3), mean, na.rm = TRUE)  # [N x T] origin margin
oa_X1 <- apply(X1, c(1, 3), mean, na.rm = TRUE)
oa_X2 <- apply(X2, c(1, 3), mean, na.rm = TRUE)

da_Y  <- apply(Y,  c(2, 3), mean, na.rm = TRUE)  # [M x T] dest margin
da_X1 <- apply(X1, c(2, 3), mean, na.rm = TRUE)
da_X2 <- apply(X2, c(2, 3), mean, na.rm = TRUE)
# apply(arr, c(1,3), ...) gives [N x T] with rows=origins, cols=times.
# Access: oa_Y[i, t]  da_Y[j, t]

# ---- C. Pair-level OLS function --------------------------------------------
pair_ols <- function(i, j) {
  # R 1-indexed: response rows (t_start+1)..T_, lag rows t_start..(T_-1)
  ts_r <- t_start + 1L
  te_r <- T_
  Te   <- te_r - ts_r + 1L

  y_reg <- Y [i, j, ts_r:te_r]
  y_lag <- Y [i, j, t_start:(te_r - 1L)]
  x1    <- X1[i, j, ts_r:te_r]
  x2    <- X2[i, j, ts_r:te_r]

  # Build augmenting block: intercept + y_lag + x1 + x2 + 9*(p_T+1) cols
  aug_list <- vector("list", 4L + 9L * (p_T + 1L))
  aug_list[[1]] <- rep(1, Te)
  aug_list[[2]] <- y_lag
  aug_list[[3]] <- x1
  aug_list[[4]] <- x2

  col_k <- 5L
  for (l in 0:p_T) {
    s <- ts_r - l; e <- te_r - l     # R 1-indexed with lag offset l
    aug_list[[col_k]]     <- ga_Y [s:e];       col_k <- col_k + 1L
    aug_list[[col_k]]     <- ga_X1[s:e];       col_k <- col_k + 1L
    aug_list[[col_k]]     <- ga_X2[s:e];       col_k <- col_k + 1L
    aug_list[[col_k]]     <- oa_Y [i, s:e];    col_k <- col_k + 1L
    aug_list[[col_k]]     <- oa_X1[i, s:e];    col_k <- col_k + 1L
    aug_list[[col_k]]     <- oa_X2[i, s:e];    col_k <- col_k + 1L
    aug_list[[col_k]]     <- da_Y [j, s:e];    col_k <- col_k + 1L
    aug_list[[col_k]]     <- da_X1[j, s:e];    col_k <- col_k + 1L
    aug_list[[col_k]]     <- da_X2[j, s:e];    col_k <- col_k + 1L
  }

  A    <- do.call(cbind, aug_list)
  # FIX: combined NaN mask — exclude rows where y_reg OR any A column is NA
  mask <- complete.cases(cbind(y_reg, A))
  if (sum(mask) < (ncol(A) + 2L)) return(NULL)

  cf <- tryCatch(
    .lm.fit(A[mask, , drop = FALSE], y_reg[mask])$coefficients,
    error = function(e) NULL
  )
  if (is.null(cf) || !all(is.finite(cf))) return(NULL)

  # cf[1]=intercept, cf[2]=rho, cf[3]=beta1, cf[4]=beta2
  list(rho = pmax(pmin(cf[2], 0.99), -0.99),
       b1  = cf[3],
       b2  = cf[4])
}

# ---- D. Run all N×M pairs --------------------------------------------------
cat("\nEstimating all pairs...\n")
b_rho  <- matrix(NA_real_, N, M)
b_b1   <- matrix(NA_real_, N, M)
b_b2   <- matrix(NA_real_, N, M)
n_valid <- 0L

for (i in seq_len(N)) {
  for (j in seq_len(M)) {
    res <- tryCatch(pair_ols(i, j), error = function(e) NULL)
    if (!is.null(res)) {
      b_rho[i, j] <- res$rho
      b_b1 [i, j] <- res$b1
      b_b2 [i, j] <- res$b2
      n_valid      <- n_valid + 1L
    }
  }
  if (i %% 5 == 0) cat(sprintf("  Origin %d/%d...\n", i, N))
}
cat(sprintf("Valid pairs: %d / %d (%.1f%%)\n", n_valid, N * M,
            100 * n_valid / (N * M)))

# ---- E. Mean group estimates -----------------------------------------------
rho_MG  <- mean(b_rho, na.rm = TRUE)
b1_MG   <- mean(b_b1,  na.rm = TRUE)
b2_MG   <- mean(b_b2,  na.rm = TRUE)
rho_c   <- pmax(pmin(b_rho, 0.99), -0.99)
theta1  <- mean(b_b1 / (1 - rho_c), na.rm = TRUE)   # long-run GDP_origin
theta2  <- mean(b_b2 / (1 - rho_c), na.rm = TRUE)   # long-run GDP_dest

cat(sprintf("\nMean group estimates:\n"))
cat(sprintf("  rho_MG  = %.4f\n", rho_MG))
cat(sprintf("  beta1_MG = %.4f\n", b1_MG))
cat(sprintf("  beta2_MG = %.4f\n", b2_MG))
cat(sprintf("  theta1_MG = %.4f\n", theta1))
cat(sprintf("  theta2_MG = %.4f\n", theta2))

# ---- F. Margin-based variance estimator (Theorem 2) -----------------------
# SE^2 = V_I + V_J = var(b_i.) / N + var(b_.j) / M
mvar <- function(b_mat) {
  b_i <- rowMeans(b_mat, na.rm = TRUE)   # [N] origin-margin means
  b_j <- colMeans(b_mat, na.rm = TRUE)   # [M] dest-margin means
  V_I <- var(b_i, na.rm = TRUE) / N
  V_J <- var(b_j, na.rm = TRUE) / M
  V_I + V_J
}

se2_rho <- mvar(b_rho);  se2_b1 <- mvar(b_b1);  se2_b2 <- mvar(b_b2)

# Delta method for long-run standard errors
# d(theta)/d(beta) = 1/(1-rho);  d(theta)/d(rho) = beta/(1-rho)^2
g1_b1  <- 1 / (1 - rho_MG);  g1_rho <- b1_MG / (1 - rho_MG)^2
g2_b2  <- 1 / (1 - rho_MG);  g2_rho <- b2_MG / (1 - rho_MG)^2

se2_th1 <- g1_b1^2 * se2_b1 + g1_rho^2 * se2_rho
se2_th2 <- g2_b2^2 * se2_b2 + g2_rho^2 * se2_rho

se_rho <- sqrt(max(se2_rho, 1e-12));  se_b1  <- sqrt(max(se2_b1, 1e-12))
se_b2  <- sqrt(max(se2_b2, 1e-12));   se_th1 <- sqrt(max(se2_th1, 1e-12))
se_th2 <- sqrt(max(se2_th2, 1e-12))

pval_fn  <- function(e, s) 2 * pnorm(-abs(e / s))
stars_fn <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**",
                         ifelse(p < 0.10, "*", "")))

# ---- G. Compile results table ----------------------------------------------
results <- data.frame(
  coef      = c("rho_MG", "beta1_MG", "beta2_MG", "theta1_MG", "theta2_MG"),
  label     = c("Persistence", "Short-run GDP, origin", "Short-run GDP, dest",
                "Long-run GDP, origin", "Long-run GDP, dest"),
  estimate  = round(c(rho_MG, b1_MG,  b2_MG,  theta1,  theta2),  4),
  se        = round(c(se_rho, se_b1,  se_b2,  se_th1,  se_th2),  4),
  tstat     = round(c(rho_MG / se_rho, b1_MG / se_b1, b2_MG / se_b2,
                       theta1 / se_th1, theta2 / se_th2), 3),
  pval      = round(pval_fn(
    c(rho_MG, b1_MG, b2_MG, theta1, theta2),
    c(se_rho, se_b1, se_b2, se_th1, se_th2)), 4),
  stringsAsFactors = FALSE
)
results$stars <- stars_fn(results$pval)
results$n_valid_pairs <- n_valid

cat("\n=== ESTIMATION RESULTS ===\n")
print(results[, c("label","estimate","se","tstat","pval","stars")], row.names = FALSE)

write.csv(results, "estimation_results.csv", row.names = FALSE)
cat("\nSaved: estimation_results.csv\n")
cat("E03_estimation.R done.  Run E04_results_excel.R next.\n")
