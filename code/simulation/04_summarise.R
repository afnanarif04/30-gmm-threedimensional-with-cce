# =============================================================================
# 04_summarise.R — Summarise Monte Carlo Results into Tables 1 and 2
# Dynamic 3D-CCE Mean Group Estimator
# Author: [anonymous]
# =============================================================================
# INPUT:  sim_results_C1.rds ... sim_results_C4.rds (from 03_simulation.R)
# OUTPUT: mc_table1_bias_rmse.csv, mc_table2_size_coverage.csv
# =============================================================================

if (interactive() && requireNamespace("rstudioapi", quietly=TRUE) &&
    rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ALPHA <- 0.05
GRID  <- data.frame(
  config = c("C1","C2","C3","C4"),
  label  = c("(20, 20, 50)","(25, 25, 100)","(40, 40, 100)","(50, 50, 150)"),
  stringsAsFactors = FALSE
)

t1_rows <- list()
t2_rows <- list()

for (k in seq_len(nrow(GRID))) {
  cfg      <- GRID[k, ]
  rds_file <- sprintf("sim_results_%s.rds", cfg$config)
  if (!file.exists(rds_file)) {
    cat(sprintf("[%s] rds not found — skipping\n", cfg$config)); next
  }

  mat   <- readRDS(rds_file)
  valid <- complete.cases(mat[, c("dyn_rho","dyn_beta","dyn_theta"), drop = FALSE])
  m     <- mat[valid, , drop = FALSE]
  n_v   <- sum(valid)
  cat(sprintf("[%s] valid=%d/%d\n", cfg$config, n_v, nrow(mat)))

  # ---- Table 1: Bias and RMSE ----------------------------------------------
  bias_rho   <- mean(m[,"dyn_rho"]   - m[,"true_rho"])
  rmse_rho   <- sqrt(mean((m[,"dyn_rho"]   - m[,"true_rho"])^2))
  bias_beta  <- mean(m[,"dyn_beta"]  - m[,"true_beta"])
  rmse_beta  <- sqrt(mean((m[,"dyn_beta"]  - m[,"true_beta"])^2))
  bias_theta <- mean(m[,"dyn_theta"] - m[,"true_theta"])
  rmse_theta <- sqrt(mean((m[,"dyn_theta"] - m[,"true_theta"])^2))

  t1_rows[[k]] <- data.frame(
    `(N, M, T)`          = cfg$label,
    `Bias (persistence)` = round(bias_rho,   4),
    `RMSE (persistence)` = round(rmse_rho,   4),
    `Bias (short run)`   = round(bias_beta,  4),
    `RMSE (short run)`   = round(rmse_beta,  4),
    `Bias (long run)`    = round(bias_theta, 4),
    `RMSE (long run)`    = round(rmse_theta, 4),
    `N valid`            = n_v,
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # ---- Table 2: Size and coverage ------------------------------------------
  # Dynamic 3D: t-stat using margin-based SE
  t_dyn  <- (m[,"dyn_rho"] - m[,"true_rho"]) / pmax(m[,"dyn_se_rho"], 1e-8)
  size_dyn <- mean(abs(t_dyn) > qnorm(1 - ALPHA/2), na.rm = TRUE)
  cov_dyn  <- mean(abs(t_dyn) <= qnorm(1 - ALPHA/2), na.rm = TRUE)

  # Static 3D: proxy size via empirical SD
  sd_sta   <- sd(m[,"sta_beta"], na.rm = TRUE)
  size_sta <- mean(abs((m[,"sta_beta"] - 1.0) / max(sd_sta, 1e-8)) >
                     qnorm(1 - ALPHA/2), na.rm = TRUE)

  # Naive 2D
  sd_nd   <- sd(m[,"nd_rho"], na.rm = TRUE)
  size_nd <- mean(abs((m[,"nd_rho"] - m[,"true_rho"]) / max(sd_nd, 1e-8)) >
                    qnorm(1 - ALPHA/2), na.rm = TRUE)

  # Pooled FE
  sd_fe   <- sd(m[,"fe_beta"], na.rm = TRUE)
  size_fe <- mean(abs((m[,"fe_beta"] - 1.0) / max(sd_fe, 1e-8)) >
                    qnorm(1 - ALPHA/2), na.rm = TRUE)

  t2_rows[[k]] <- data.frame(
    `(N, M, T)`              = cfg$label,
    `Dynamic 3D size`        = round(size_dyn, 3),
    `Dynamic 3D coverage`    = round(cov_dyn,  3),
    `Static 3D size`         = round(size_sta, 3),
    `Naive 2D dynamic size`  = round(size_nd,  3),
    `Pooled FE size`         = round(size_fe,  3),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

table1 <- do.call(rbind, t1_rows)
table2 <- do.call(rbind, t2_rows)

write.csv(table1, "mc_table1_bias_rmse.csv",    row.names = FALSE)
write.csv(table2, "mc_table2_size_coverage.csv", row.names = FALSE)
cat("Saved mc_table1_bias_rmse.csv and mc_table2_size_coverage.csv\n")
cat("04_summarise.R done. Run 05_empirical_data.R next.\n")
