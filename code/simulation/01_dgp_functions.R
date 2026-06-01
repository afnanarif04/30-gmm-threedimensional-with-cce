# =============================================================================
# 01_dgp_functions.R — Data Generating Process
# Dynamic 3D-CCE Mean Group Estimator
# Paper: [anonymous]
# Author: [anonymous]
# =============================================================================
# Hierarchical DGP:
#   y_{ijt} = rho_{ij} * y_{ij,t-1} + beta_{ij} * x_{ijt} + u_{ijt}
#   u_{ijt} = lambda_G' f_t + lambda_I' g_{it} + lambda_J' h_{jt} + eps_{ijt}
# Slopes: beta_{ij} = beta0 + nu_I_i + nu_J_j + nu_P_ij (hierarchical)
#         rho_{ij}  = rho0  + nu_rI_i + nu_rJ_j (clipped to [0.1, 0.8])
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly=TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

generate_3d_panel <- function(N, M, T,
    beta0  = 1.0,
    rho0   = 0.4,
    sig_I  = 0.15,   # SD of origin-margin slope component
    sig_J  = 0.15,   # SD of destination-margin slope component
    sig_P  = 0.05,   # SD of pair-specific slope component
    sig_rI = 0.05,   # SD of origin-margin rho component
    sig_rJ = 0.05,   # SD of destination-margin rho component
    mG     = 2,      # number of global factors
    mI     = 2,      # number of origin-specific factors
    mJ     = 2,      # number of destination-specific factors
    ar_f   = 0.5,    # AR(1) coefficient of all factor processes
    T_burn = 30      # burn-in periods
) {

  T_total <- T + T_burn

  # ---- Heterogeneous slopes (hierarchical decomposition) -------------------
  nu_I  <- rnorm(N, 0, sig_I)
  nu_J  <- rnorm(M, 0, sig_J)
  nu_P  <- matrix(rnorm(N * M, 0, sig_P), nrow = N, ncol = M)

  beta_ij <- matrix(beta0, N, M) +
             nu_I %*% matrix(1, 1, M) +
             matrix(1, N, 1) %*% t(nu_J) + nu_P

  nu_rI <- rnorm(N, 0, sig_rI)
  nu_rJ <- rnorm(M, 0, sig_rJ)
  rho_ij <- pmax(pmin(
    matrix(rho0, N, M) +
    nu_rI %*% matrix(1, 1, M) +
    matrix(1, N, 1) %*% matrix(nu_rJ, nrow = 1),
    0.99), 0.01)

  # ---- AR(1) factor simulation ---------------------------------------------
  sim_ar <- function(n_fac) {
    F  <- matrix(0, T_total, n_fac)
    F[1, ] <- rnorm(n_fac)
    sd_inn <- sqrt(1 - ar_f^2)
    for (t in 2:T_total)
      F[t, ] <- ar_f * F[t-1, ] + rnorm(n_fac, 0, sd_inn)
    F
  }

  F_G <- sim_ar(mG)                              # global factors  [T_total x mG]
  G_O <- lapply(seq_len(N), function(i) sim_ar(mI))  # origin factors
  H_D <- lapply(seq_len(M), function(j) sim_ar(mJ))  # dest factors

  # ---- Loadings (positive mean ensures rank condition A3) ------------------
  lG <- array(runif(N*M*mG, 0.5, 1.5), c(N, M, mG))
  gG <- array(runif(N*M*mG, 0.5, 1.5), c(N, M, mG))
  lI <- array(runif(N*M*mI, 0.5, 1.5), c(N, M, mI))
  gI <- array(runif(N*M*mI, 0.5, 1.5), c(N, M, mI))
  lJ <- array(runif(N*M*mJ, 0.5, 1.5), c(N, M, mJ))
  gJ <- array(runif(N*M*mJ, 0.5, 1.5), c(N, M, mJ))

  # ---- Simulate Y and X arrays [N x M x T_total] --------------------------
  Y <- array(0, c(N, M, T_total))
  X <- array(0, c(N, M, T_total))

  for (t in seq_len(T_total)) {
    for (i in seq_len(N)) {
      for (j in seq_len(M)) {
        X[i, j, t] <- sum(gG[i,j,] * F_G[t,]) +
                      sum(gI[i,j,] * G_O[[i]][t,]) +
                      sum(gJ[i,j,] * H_D[[j]][t,]) + rnorm(1, 0, 0.5)
      }
    }
  }

  Y[,, 1] <- matrix(rnorm(N * M), N, M)
  for (t in 2:T_total) {
    for (i in seq_len(N)) {
      for (j in seq_len(M)) {
        u <- sum(lG[i,j,] * F_G[t,]) +
             sum(lI[i,j,] * G_O[[i]][t,]) +
             sum(lJ[i,j,] * H_D[[j]][t,]) + rnorm(1, 0, 0.5)
        Y[i, j, t] <- rho_ij[i,j] * Y[i,j,t-1] + beta_ij[i,j] * X[i,j,t] + u
      }
    }
  }

  # ---- Drop burn-in ---------------------------------------------------------
  ti <- (T_burn + 1):(T_burn + T)
  list(
    Y       = Y[,, ti],
    X       = X[,, ti],
    beta_ij = beta_ij,
    rho_ij  = rho_ij,
    beta0   = beta0,
    rho0    = mean(rho_ij),
    theta0  = mean(beta_ij / (1 - rho_ij))
  )
}

cat("01_dgp_functions.R loaded. Run 02_estimators.R next.\n")
