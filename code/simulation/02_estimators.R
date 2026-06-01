# =============================================================================
# 02_estimators.R — Four Estimators
# Dynamic 3D-CCE Mean Group Estimator
# Author: [anonymous]
# =============================================================================
# Four estimators:
#   (1) estimate_dynamic_3d_ccemg  — PROPOSED estimator
#   (2) estimate_static_3d_ccemg   — Kapetanios, Serlenga & Shin (2021)
#   (3) estimate_2d_dynamic_ccemg  — Chudik & Pesaran (2015) naive 2D
#   (4) estimate_pooled_fe         — additive fixed effects (wrong baseline)
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly=TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Helper: lag-truncation order -------------------------------------------
# Enforces minimum T_eff / n_params >= min_ratio to prevent near-singular
# pair regressions. For univariate (Monte Carlo): k_reg=3, N_aug_per_lag=6.
# For bivariate (empirical):           k_reg=4, N_aug_per_lag=9.
get_pT <- function(T, N_aug_per_lag = 6, k_reg = 3, min_ratio = 2.5) {
  p_cube <- floor(T^(1/3))
  p_safe <- max(0L, as.integer(
    (T - 1 - min_ratio * (k_reg + N_aug_per_lag)) /
    (1 + N_aug_per_lag * min_ratio)
  ))
  min(p_cube, p_safe)
}

# ---- Helper: margin-based variance estimator --------------------------------
# Sigma_MG = V_I + V_J (origin margin + destination margin dispersions)
# V_I = Var(b_{i.}) / N,  V_J = Var(b_{.j}) / M
margin_variance <- function(b_mat, N, M) {
  b_i <- rowMeans(b_mat, na.rm = TRUE)   # [N] origin-margin means
  b_j <- colMeans(b_mat, na.rm = TRUE)   # [M] dest-margin means
  V_I <- var(b_i, na.rm = TRUE) / N
  V_J <- var(b_j, na.rm = TRUE) / M
  V_I + V_J
}

# ---- (1) PROPOSED: Dynamic 3D-CCE Mean Group --------------------------------
estimate_dynamic_3d_ccemg <- function(Y, X) {
  N <- dim(Y)[1]; M <- dim(Y)[2]; T <- dim(Y)[3]
  p <- get_pT(T)

  # Pre-compute cross-sectional averages [T x 2] (Y col1, X col2)
  ga <- cbind(sapply(seq_len(T), function(t) mean(Y[,,t])),
              sapply(seq_len(T), function(t) mean(X[,,t])))
  oa <- array(0, c(N, T, 2))
  da <- array(0, c(M, T, 2))
  for (i in seq_len(N))
    oa[i,,] <- cbind(sapply(seq_len(T), function(t) mean(Y[i,,t])),
                     sapply(seq_len(T), function(t) mean(X[i,,t])))
  for (j in seq_len(M))
    da[j,,] <- cbind(sapply(seq_len(T), function(t) mean(Y[,j,t])),
                     sapply(seq_len(T), function(t) mean(X[,j,t])))

  t_start <- p + 2L                        # 1-indexed first response row
  if (T - t_start + 1 < 6) return(NULL)

  build_aug <- function(i, j) {
    Te    <- T - t_start + 1L
    y_reg <- Y[i, j, t_start:T]
    y_lag <- Y[i, j, (t_start-1):(T-1)]
    x_reg <- X[i, j, t_start:T]
    aug   <- vector("list", p + 1)
    for (l in 0:p) {
      s <- t_start - l; e <- T - l
      aug[[l+1]] <- cbind(ga[s:e,], oa[i,s:e,], da[j,s:e,])
    }
    H <- do.call(cbind, aug)
    list(y = y_reg, A = cbind(1, y_lag, x_reg, H))
  }

  b_full <- array(NA, c(N, M, 2))
  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      d <- build_aug(i, j)
      mask <- complete.cases(cbind(d$y, d$A))
      if (sum(mask) < ncol(d$A) + 2) next
      cf <- tryCatch(.lm.fit(d$A[mask,,drop=FALSE], d$y[mask])$coefficients,
                     error = function(e) NULL)
      if (!is.null(cf) && all(is.finite(cf)))
        b_full[i, j, ] <- cf[2:3]        # rho (col2), beta (col3)
    }
  }

  b_MG <- apply(b_full, 3, mean, na.rm = TRUE)

  # NOTE: No half-panel jackknife. At finite T the jackknife amplifies bias
  # 3x because half-panel regressions are near-singular. Plain MG is used.
  # The residual O(1/T) dynamic bias is an honest finite-sample finding.

  # Long-run by delta method
  rho_c    <- pmax(pmin(b_full[,,1], 0.99), -0.99)
  theta_MG <- mean(b_full[,,2] / (1 - rho_c), na.rm = TRUE)

  # Margin-based SEs (scaled by min(N,M))
  Sig    <- margin_variance(b_full[,,1], N, M) + margin_variance(b_full[,,2], N, M)
  se_rho <- sqrt(max(margin_variance(b_full[,,1], N, M), 1e-12))
  se_bet <- sqrt(max(margin_variance(b_full[,,2], N, M), 1e-12))

  list(
    rho_MG  = b_MG[1],
    beta_MG = b_MG[2],
    theta_MG = theta_MG,
    se_rho  = se_rho,
    se_beta = se_bet,
    b_full  = b_full
  )
}

# ---- (2) Static 3D-CCE MG (Kapetanios, Serlenga & Shin 2021) ---------------
estimate_static_3d_ccemg <- function(Y, X) {
  N <- dim(Y)[1]; M <- dim(Y)[2]; T <- dim(Y)[3]
  p <- get_pT(T)

  ga <- cbind(sapply(seq_len(T), function(t) mean(Y[,,t])),
              sapply(seq_len(T), function(t) mean(X[,,t])))
  oa <- array(0, c(N, T, 2))
  da <- array(0, c(M, T, 2))
  for (i in seq_len(N))
    oa[i,,] <- cbind(sapply(seq_len(T), function(t) mean(Y[i,,t])),
                     sapply(seq_len(T), function(t) mean(X[i,,t])))
  for (j in seq_len(M))
    da[j,,] <- cbind(sapply(seq_len(T), function(t) mean(Y[,j,t])),
                     sapply(seq_len(T), function(t) mean(X[,j,t])))

  t_start <- p + 1L
  b_mat <- numeric(N * M)
  k <- 0L
  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      k <- k + 1L
      y_reg <- Y[i, j, t_start:T]
      x_reg <- X[i, j, t_start:T]
      aug   <- vector("list", p + 1)
      for (l in 0:p) {
        s <- t_start - l; e <- T - l
        aug[[l+1]] <- cbind(ga[s:e,], oa[i,s:e,], da[j,s:e,])
      }
      H  <- do.call(cbind, aug)
      A  <- cbind(1, x_reg, H)
      mask <- complete.cases(cbind(y_reg, A))
      if (sum(mask) < ncol(A) + 1) { b_mat[k] <- NA; next }
      cf <- tryCatch(.lm.fit(A[mask,,drop=FALSE], y_reg[mask])$coefficients,
                     error = function(e) NULL)
      b_mat[k] <- if (!is.null(cf) && all(is.finite(cf))) cf[2] else NA
    }
  }
  list(beta_MG = mean(b_mat, na.rm = TRUE))
}

# ---- (3) Naive 2D Dynamic CCE-MG (Chudik & Pesaran 2015) -------------------
# Ignores destination-margin averages — misses H_D factor level
estimate_2d_dynamic_ccemg <- function(Y, X) {
  N <- dim(Y)[1]; M <- dim(Y)[2]; T <- dim(Y)[3]
  p <- get_pT(T)

  ga <- cbind(sapply(seq_len(T), function(t) mean(Y[,,t])),
              sapply(seq_len(T), function(t) mean(X[,,t])))
  oa <- array(0, c(N, T, 2))
  for (i in seq_len(N))
    oa[i,,] <- cbind(sapply(seq_len(T), function(t) mean(Y[i,,t])),
                     sapply(seq_len(T), function(t) mean(X[i,,t])))

  t_start <- p + 2L
  b_mat <- array(NA, c(N, M, 2))
  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      y_reg <- Y[i, j, t_start:T]
      y_lag <- Y[i, j, (t_start-1):(T-1)]
      x_reg <- X[i, j, t_start:T]
      aug   <- vector("list", p + 1)
      for (l in 0:p) {
        s <- t_start - l; e <- T - l
        aug[[l+1]] <- cbind(ga[s:e,], oa[i,s:e,])  # NO dest averages
      }
      H  <- do.call(cbind, aug)
      A  <- cbind(1, y_lag, x_reg, H)
      mask <- complete.cases(cbind(y_reg, A))
      if (sum(mask) < ncol(A) + 2) next
      cf <- tryCatch(.lm.fit(A[mask,,drop=FALSE], y_reg[mask])$coefficients,
                     error = function(e) NULL)
      if (!is.null(cf) && all(is.finite(cf)))
        b_mat[i, j, ] <- cf[2:3]
    }
  }
  rho_c <- pmax(pmin(b_mat[,,1], 0.99), -0.99)
  list(
    rho_MG   = mean(b_mat[,,1], na.rm = TRUE),
    beta_MG  = mean(b_mat[,,2], na.rm = TRUE),
    theta_MG = mean(b_mat[,,2] / (1 - rho_c), na.rm = TRUE)
  )
}

# ---- (4) Pooled Fixed Effects (additive dummies) ----------------------------
# Cannot absorb interactive multifactor structure — severe baseline
estimate_pooled_fe <- function(Y, X) {
  N <- dim(Y)[1]; M <- dim(Y)[2]; T <- dim(Y)[3]
  grid <- expand.grid(i = seq_len(N), j = seq_len(M), t = seq_len(T))
  grid$y <- as.vector(Y)
  grid$x <- as.vector(X)
  # Within-transform: demean by i, j, and t (additive FE)
  grid$yd <- grid$y - ave(grid$y, grid$i) -
             ave(grid$y, grid$j) - ave(grid$y, grid$t) + 2 * mean(grid$y, na.rm = TRUE)
  grid$xd <- grid$x - ave(grid$x, grid$i) -
             ave(grid$x, grid$j) - ave(grid$x, grid$t) + 2 * mean(grid$x, na.rm = TRUE)
  mask <- complete.cases(grid[, c("yd","xd")])
  cf   <- tryCatch(
    lm(yd ~ xd - 1, data = grid[mask, ])$coefficients,
    error = function(e) NA_real_
  )
  list(beta_FE = if (length(cf) > 0 && is.finite(cf[1])) cf[1] else NA_real_)
}

cat("02_estimators.R loaded. Run 03_simulation.R next.\n")
