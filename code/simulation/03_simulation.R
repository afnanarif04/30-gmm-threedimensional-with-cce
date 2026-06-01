# =============================================================================
# 03_simulation.R — Main Monte Carlo Simulation
# Dynamic 3D-CCE Mean Group Estimator
# Author: [anonymous]
# =============================================================================
# GRID: 4 cells x 2000 replications
#   C1: (N=20, M=20, T=50)   C2: (N=25, M=25, T=100)
#   C3: (N=40, M=40, T=100)  C4: (N=50, M=50, T=150)
# OUTPUT: sim_results_C1.rds ... sim_results_C4.rds
#         mc_table1_bias_rmse.csv, mc_table2_size_coverage.csv
# RUNTIME: ~1h with 15 cores, ~8h single-core
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly=TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(parallel, quietly = TRUE)
source("01_dgp_functions.R")
source("02_estimators.R")

# ---- Parameters ------------------------------------------------------------
R_REPS  <- 2000
N_CORES <- max(1L, parallel::detectCores() - 1L)
BETA0   <- 1.0
RHO0    <- 0.4
ALPHA   <- 0.05
SEED    <- 20250529L

GRID <- data.frame(
  config = c("C1","C2","C3","C4"),
  N      = c(20L, 25L, 40L, 50L),
  M      = c(20L, 25L, 40L, 50L),
  T      = c(50L, 100L, 100L, 150L),
  label  = c("(20,20,50)","(25,25,100)","(40,40,100)","(50,50,150)"),
  stringsAsFactors = FALSE
)

RESULT_COLS <- c(
  "true_rho","true_beta","true_theta",
  "dyn_rho","dyn_beta","dyn_theta","dyn_se_rho","dyn_reject",
  "sta_beta",
  "nd_rho","nd_beta","nd_theta",
  "fe_beta"
)

# ---- Single replication function -------------------------------------------
run_one_rep <- function(seed_val, N, M, T) {
  set.seed(seed_val)
  dat <- tryCatch(
    generate_3d_panel(N, M, T, beta0 = BETA0, rho0 = RHO0),
    error = function(e) NULL
  )
  if (is.null(dat)) return(rep(NA_real_, length(RESULT_COLS)))

  rho_true   <- mean(dat$rho_ij)
  theta_true <- dat$theta0

  num <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != 1 || !is.finite(v)) NA_real_ else v
  }

  # Proposed: dynamic 3D-CCE MG
  e1 <- tryCatch(estimate_dynamic_3d_ccemg(dat$Y, dat$X), error = function(e) NULL)
  dyn_reject <- if (!is.null(e1) && !is.null(e1$se_rho) && is.finite(e1$se_rho)) {
    as.numeric(abs((e1$rho_MG - rho_true) / max(e1$se_rho, 1e-8)) > qnorm(1 - ALPHA/2))
  } else NA_real_

  # Static 3D-CCE MG
  e2 <- tryCatch(estimate_static_3d_ccemg(dat$Y, dat$X),    error = function(e) NULL)
  # Naive 2D dynamic CCE-MG
  e3 <- tryCatch(estimate_2d_dynamic_ccemg(dat$Y, dat$X),   error = function(e) NULL)
  # Pooled FE
  e4 <- tryCatch(estimate_pooled_fe(dat$Y, dat$X),          error = function(e) NULL)

  g  <- function(e, nm) if (!is.null(e) && !is.null(e[[nm]])) e[[nm]] else NA_real_

  c(
    true_rho   = num(rho_true),
    true_beta  = num(BETA0),
    true_theta = num(theta_true),
    dyn_rho    = num(g(e1, "rho_MG")),
    dyn_beta   = num(g(e1, "beta_MG")),
    dyn_theta  = num(g(e1, "theta_MG")),
    dyn_se_rho = num(g(e1, "se_rho")),
    dyn_reject = num(dyn_reject),
    sta_beta   = num(g(e2, "beta_MG")),
    nd_rho     = num(g(e3, "rho_MG")),
    nd_beta    = num(g(e3, "beta_MG")),
    nd_theta   = num(g(e3, "theta_MG")),
    fe_beta    = num(g(e4, "beta_FE"))
  )
}

# ---- Main loop --------------------------------------------------------------
set.seed(SEED)
cat(sprintf("Monte Carlo: %d reps x %d cells x %d cores\n",
            R_REPS, nrow(GRID), N_CORES))

for (k in seq_len(nrow(GRID))) {
  cfg   <- GRID[k, ]
  N <- cfg$N; M <- cfg$M; T <- cfg$T; label <- cfg$label
  rds_file <- sprintf("sim_results_%s.rds", cfg$config)

  if (file.exists(rds_file)) {
    cat(sprintf("[%s] %s — cache found, skipping\n", cfg$config, label))
    next
  }

  cat(sprintf("[%s] %s — starting %d reps ...\n", cfg$config, label, R_REPS))
  t0 <- proc.time()[3]
  rep_seeds <- sample.int(.Machine$integer.max, R_REPS)

  if (N_CORES > 1L) {
    cl <- makeCluster(N_CORES)
    proj_dir <- getwd()
    clusterExport(cl, "proj_dir", envir = environment())
    invisible(clusterEvalQ(cl, {
      setwd(proj_dir)
      suppressMessages({
        source("01_dgp_functions.R")
        source("02_estimators.R")
      })
      TRUE
    }))
    clusterExport(cl, c("run_one_rep","RESULT_COLS","BETA0","RHO0","ALPHA"),
                  envir = .GlobalEnv)
    clusterExport(cl, c("N","M","T"), envir = environment())
    raw_list <- parLapply(cl, rep_seeds, function(s) run_one_rep(s, N, M, T))
    stopCluster(cl)
  } else {
    raw_list <- lapply(seq_along(rep_seeds), function(r) {
      if (r %% 200 == 0) cat(sprintf("  rep %d/%d\n", r, R_REPS))
      run_one_rep(rep_seeds[r], N, M, T)
    })
  }

  results_mat <- matrix(unlist(raw_list), nrow = R_REPS,
                        ncol = length(RESULT_COLS), byrow = TRUE)
  colnames(results_mat) <- RESULT_COLS
  saveRDS(results_mat, rds_file)
  cat(sprintf("  Done in %.0fs — saved %s\n", proc.time()[3] - t0, rds_file))
}

cat("03_simulation.R complete. Run 04_summarise.R next.\n")
