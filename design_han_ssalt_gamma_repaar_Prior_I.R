rm(list = ls())

library(parallel)
library(doParallel)
library(foreach)
library(cmdstanr)
library(ggplot2)
library(cowplot)
library(lattice)

par(cex = 1.25, mar = c(4.75, 4.85, 2.25, 1.5))

source("/mnt/nfs/home/nkp117/largefiles/comet_test/all_functions.R")

# ============================================================
# Fixed stress boundaries
# ============================================================

s0 <- 1/293    # normal use stress  (fixed lower bound)
s2 <- 1/353    # highest accelerated stress (fixed upper bound)

# ============================================================
# s1 grid: x1 = (s1-s0)/(s2-s0) from 0.1 to 0.7
# Avoid 0 and 1 — those collapse the stress contrast
# ============================================================

x1_grid  <- seq(0.1, 0.9, length.out = 9)   
s1_grid  <- s0 + x1_grid * (s2 - s0)

cat("s1 grid (1/s1 rounded):\n"); print(round(1/s1_grid))
cat("\nx1 grid:\n");               print(round(x1_grid, 3))
cat("\n")

# ============================================================
# Master output directory
# ============================================================

master_dir <- "/mnt/nfs/home/nkp117/largefiles/Prior_I/"
dir.create(master_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Fixed model parameters
# ============================================================

set.seed(2026)

n  <- 35
tc <- 6
p  <- 0.1
x0 <- 0
x2 <- 1

para_hat  <- c(4.5064079, -4.7131110, 0.7692292,
               2.0409831, -1.2277314, 1.5320872)
para_true <- para_hat

# ============================================================
# Gamma hyperparameters
# ============================================================

gamma_hyper <- function(mean_sd) {
  mean_x <- mean_sd[1];  sd_x <- mean_sd[2]
  c(mean_x^2 / sd_x^2,  mean_x / sd_x^2)
}

hpar11 <- c(0.1634133, 0.3704860);  hpar21 <- c(4.280502, 1.2736822);  hpar31 <- c(1.200586, 1.2724302)
hpar12 <- c(0.1526861, 0.1549505);  hpar22 <- c(1.402464, 0.5039255);  hpar32 <- c(1.698900, 0.4604109)

gphi11 <- gamma_hyper(hpar11);  gphi21 <- gamma_hyper(hpar21);  gphi31 <- gamma_hyper(hpar31)
gphi12 <- gamma_hyper(hpar12);  gphi22 <- gamma_hyper(hpar22);  gphi32 <- gamma_hyper(hpar32)

init_fun <- function() {
  list(
    phi_11 = max(rgamma(1, gphi11[1], gphi11[2]), 1e-2),
    phi_21 = max(rgamma(1, gphi21[1], gphi21[2]), 1e-2),
    phi_31 = max(rgamma(1, gphi31[1], gphi31[2]), 1e-2),
    phi_12 = max(rgamma(1, gphi12[1], gphi12[2]), 1e-2),
    phi_22 = max(rgamma(1, gphi22[1], gphi22[2]), 1e-2),
    phi_32 = max(rgamma(1, gphi32[1], gphi32[2]), 1e-2)
  )
}

# ============================================================
# Diagnostics thresholds
# ============================================================

SD_THRESHOLD   <- 1e4
RHAT_THRESHOLD <- 1.01
ESS_THRESHOLD  <- 300
WARN_VALID     <- 950
ERROR_VALID    <- 900

# ============================================================
# Design grid and B
# ============================================================

#tau_grid <- seq(1.59375, 5.75, length.out = 20)      # 25 tau values
#tau_grid <- seq(0.05, 5.95, length.out = 25)
tau_grid <- seq(0.05, 5.95, length.out = 25)
B        <- 1000

# ============================================================
# Parallelisation — setup ONCE, reused for all loops
# ============================================================

total_cores    <- as.integer(Sys.getenv("SLURM_NTASKS", unset = detectCores()))
chains_per_fit <- 3
n_workers      <- max(1, floor(total_cores / chains_per_fit))

cat("Detected cores:", total_cores, "\n")
cat("Chains per fit:", chains_per_fit, "\n")
cat("Using workers:",  n_workers, "\n\n")

cl <- makePSOCKcluster(n_workers)
registerDoParallel(cl)

# ============================================================
# Compile Stan model ONCE — shared across all s1 and tau
# ============================================================

stan_file <- file.path(cmdstan_path(),
                       "examples/bernoulli/bernoulli_han_ssalt_gamma_repara.stan")
mod <- cmdstan_model(stan_file, compile = TRUE)

# ============================================================
# Storage for combined surface — accumulated across s1 loop
# ============================================================

surface_summary   <- data.frame()
surface_avg_sds   <- data.frame()
surface_avg_means <- data.frame()

# ============================================================
# OUTER LOOP: over s1 values
# ============================================================

for (s1_idx in seq_along(s1_grid)) {

  s1       <- s1_grid[s1_idx]
  x1       <- x1_grid[s1_idx]
  s1_label <- as.character(round(1/s1))

  out_dir <- paste0(master_dir, "s1_", s1_label, "/")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n################################################################\n")
  cat(sprintf("## s1 = %d / %d :  1/s1 = %s  (x1 = %.3f)\n",
              s1_idx, length(s1_grid), s1_label, x1))
  cat("################################################################\n\n")

  all_results <- vector("list", length(tau_grid))

  # ---- MIDDLE LOOP: tau (sequential) ----------------------

  for (i in seq_along(tau_grid)) {

    tau_val <- tau_grid[i]
    cat(sprintf("  tau = %.4f  (%d / %d)\n", tau_val, i, length(tau_grid)))

    # ---- INNER LOOP: B replications (parallel) ------------

    rep_results <- foreach(
      b              = 1:B,
      .packages      = c("cmdstanr"),
      .errorhandling = "pass"
    ) %dopar% {

      sim <- simSSAT_CR(n = n, tau = tau_val, tc = tc,
                        para_vec = para_true, x1 = x1, x2 = x2)

      t   <- pmin(sim$t, tc)
      C_i <- sim$C_i
      n1  <- sum(t <= tau_val)
      nc  <- sum(t <= tc)

      data_list <- list(
        n = n, tau = tau_val, tc = tc, n1 = n1, nc = nc,
        t = t, C_i = C_i, x1 = x1, x2 = x2, x0 = x0, p = p,
        gphi11 = gphi11, gphi21 = gphi21, gphi31 = gphi31,
        gphi12 = gphi12, gphi22 = gphi22, gphi32 = gphi32
      )

      json_file <- tempfile(fileext = ".json")
      write_stan_json(data_list, json_file)

      # --- First attempt ---
      fit <- mod$sample(
        data            = json_file,
        chains          = chains_per_fit,
        parallel_chains = chains_per_fit,
        init            = init_fun,
        iter_warmup     = 1000,
        iter_sampling   = 1000,
        seed            = 100000 + 10000*s1_idx + 100*i + b,
        refresh         = 0
      )

      diag  <- fit$diagnostic_summary()
      n_div <- sum(diag$num_divergent)
      summ  <- fit$summary()

      bad_rhat <- any(summ$rhat     > RHAT_THRESHOLD, na.rm = TRUE)
      bad_inf  <- any(is.infinite(summ$sd),            na.rm = TRUE)
      bad_sd   <- any(summ$sd       > SD_THRESHOLD,    na.rm = TRUE)
      bad_ess  <- any(summ$ess_bulk < ESS_THRESHOLD,   na.rm = TRUE) |
                  any(summ$ess_tail < ESS_THRESHOLD,   na.rm = TRUE)
      bad_fit  <- (n_div > 0) || bad_rhat || bad_inf || bad_sd || bad_ess

      # --- Refit ---
      if (bad_fit) {
        fit <- mod$sample(
          data            = json_file,
          chains          = chains_per_fit,
          parallel_chains = chains_per_fit,
          init            = init_fun,
          iter_warmup     = 2000,
          iter_sampling   = 2000,
          adapt_delta     = 0.99,
          max_treedepth   = 15,
          seed            = 200000 + 10000*s1_idx + 100*i + b,
          refresh         = 0
        )

        diag  <- fit$diagnostic_summary()
        n_div <- sum(diag$num_divergent)
        summ  <- fit$summary()

        bad_rhat <- any(summ$rhat     > RHAT_THRESHOLD, na.rm = TRUE)
        bad_inf  <- any(is.infinite(summ$sd),            na.rm = TRUE)
        bad_sd   <- any(summ$sd       > SD_THRESHOLD,    na.rm = TRUE)
        bad_ess  <- any(summ$ess_bulk < ESS_THRESHOLD,   na.rm = TRUE) |
                    any(summ$ess_tail < ESS_THRESHOLD,   na.rm = TRUE)
        bad_fit  <- (n_div > 0) || bad_rhat || bad_inf || bad_sd || bad_ess
      }

      if (bad_fit) {
        return(list(status   = "discarded", rep = b,
                    n_div    = n_div,
                    bad_rhat = bad_rhat,
                    bad_sd   = bad_sd))
      }

      list(status = "valid", rep = b,
           means  = setNames(summ$mean, summ$variable),
           sds    = setNames(summ$sd,   summ$variable),
           vars   = summ$variable)

    } # end foreach B

    valid_reps <- Filter(function(r) !inherits(r,"error") && r$status=="valid",     rep_results)
    disc_reps  <- Filter(function(r) !inherits(r,"error") && r$status=="discarded", rep_results)
    err_reps   <- Filter(function(r)  inherits(r,"error"),                           rep_results)

    n_valid <- length(valid_reps)
    cat(sprintf("    Valid: %d  Discarded: %d  Errors: %d\n",
                n_valid, length(disc_reps), length(err_reps)))

    if (n_valid < ERROR_VALID)
      stop(sprintf("FATAL: %d valid reps at s1=1/%s tau=%.4f",
                   n_valid, s1_label, tau_val))
    if (n_valid < WARN_VALID)
      warning(sprintf("Low reps: %d/%d at s1=1/%s tau=%.4f",
                      n_valid, B, s1_label, tau_val))

    mean_matrix <- do.call(rbind, lapply(valid_reps, `[[`, "means"))
    sd_matrix   <- do.call(rbind, lapply(valid_reps, `[[`, "sds"))
    colnames(mean_matrix) <- valid_reps[[1]]$vars
    colnames(sd_matrix)   <- valid_reps[[1]]$vars

    all_results[[i]] <- list(
      tau         = tau_val,
      C1_tp       = mean(sd_matrix[, "tp"]^2,     na.rm = TRUE),
      C1_logtp    = mean(sd_matrix[, "log_tp"]^2, na.rm = TRUE),
      avg_means   = colMeans(mean_matrix, na.rm = TRUE),
      avg_sds     = colMeans(sd_matrix,   na.rm = TRUE),
      n_valid     = n_valid,
      n_discarded = length(disc_reps)
    )

  } # end tau loop

  # ----------------------------------------------------------
  # Per-s1: design summary
  # ----------------------------------------------------------

  design_summary <- data.frame(
    s1          = s1,
    x1          = x1,
    tau         = sapply(all_results, `[[`, "tau"),
    C1_tp       = sapply(all_results, `[[`, "C1_tp"),
    C1_logtp    = sapply(all_results, `[[`, "C1_logtp"),
    n_valid     = sapply(all_results, `[[`, "n_valid"),
    n_discarded = sapply(all_results, `[[`, "n_discarded")
  )

  cat(sprintf("\n  Design summary s1 = 1/%s:\n", s1_label))
  print(design_summary)

  cat(sprintf("\n  Survival s1 = 1/%s:\n", s1_label))
  for (r in all_results)
    cat(sprintf("    tau=%.4f : %d/%d (%.1f%%)\n",
                r$tau, r$n_valid, B, 100*r$n_valid/B))

  # ----------------------------------------------------------
  # Per-s1: save RData
  # ----------------------------------------------------------

  save.image(file = paste0(out_dir, "preposterior_s1_", s1_label, "_full.RData"))

  # ----------------------------------------------------------
  # Per-s1: CSVs
  # ----------------------------------------------------------

  write.csv(design_summary,
            paste0(out_dir, "design_summary_s1_", s1_label, ".csv"),
            row.names = FALSE)

  valid_results <- Filter(function(r) !is.null(r$avg_sds) && length(r$avg_sds) > 1,
                          all_results)
  param_names <- names(valid_results[[1]]$avg_means)

  avg_means_df <- do.call(rbind, lapply(valid_results, function(r) {
    df <- as.data.frame(t(r$avg_means))
    df$tau <- r$tau;  df$s1 <- s1;  df$x1 <- x1;  df
  }))
  avg_sds_df <- do.call(rbind, lapply(valid_results, function(r) {
    df <- as.data.frame(t(r$avg_sds))
    df$tau <- r$tau;  df$s1 <- s1;  df$x1 <- x1;  df
  }))

  id_cols      <- c("tau","s1","x1")
  avg_means_df <- avg_means_df[, c(id_cols, setdiff(names(avg_means_df), id_cols))]
  avg_sds_df   <- avg_sds_df[,   c(id_cols, setdiff(names(avg_sds_df),   id_cols))]

  write.csv(avg_means_df,
            paste0(out_dir, "avg_means_s1_", s1_label, ".csv"), row.names = FALSE)
  write.csv(avg_sds_df,
            paste0(out_dir, "avg_sds_s1_",   s1_label, ".csv"), row.names = FALSE)

  cat(sprintf("  Saved all outputs for s1 = 1/%s  ->  %s\n", s1_label, out_dir))

  # ----------------------------------------------------------
  # Per-s1: SD and Means plots
  # ----------------------------------------------------------

  selected_params <- c("a1","a2","b1","b2","beta1","beta2",
                       "thetaa[1,1]","thetaa[2,1]",
                       "thetaa[1,2]","thetaa[2,2]","tp","log_tp")

  se_data <- do.call(rbind, lapply(valid_results, function(res)
    data.frame(tau=res$tau, parameter=param_names,
               SE=as.numeric(res$avg_sds), stringsAsFactors=FALSE)))
  mean_data <- do.call(rbind, lapply(valid_results, function(res)
    data.frame(tau=res$tau, parameter=names(res$avg_means),
               Mean=as.numeric(res$avg_means), stringsAsFactors=FALSE)))

  p_sd <- ggplot(subset(se_data, parameter %in% selected_params),
                 aes(x=tau, y=SE)) +
    geom_line(linewidth=1) + geom_point(size=2) +
    facet_wrap(~parameter, scales="free_y", ncol=4) +
    labs(x=expression(tau), y="Posterior SD",
         title=bquote("Posterior SD vs "*tau*" | s1 = 1/"*.(s1_label))) +
    theme_bw(base_size=14) + theme(plot.title=element_text(hjust=0.5))

  ggsave(paste0(out_dir,"plot_SD_s1_",   s1_label,".pdf"), p_sd,  width=14, height=10)
  ggsave(paste0(out_dir,"plot_SD_s1_",   s1_label,".eps"), p_sd,  width=14, height=10, device="eps")

  p_mean <- ggplot(subset(mean_data, parameter %in% selected_params),
                   aes(x=tau, y=Mean)) +
    geom_line(linewidth=1) + geom_point(size=2) +
    facet_wrap(~parameter, scales="free_y", ncol=4) +
    labs(x=expression(tau), y="Posterior Mean",
         title=bquote("Posterior Mean vs "*tau*" | s1 = 1/"*.(s1_label))) +
    theme_bw(base_size=14) + theme(plot.title=element_text(hjust=0.5))

  ggsave(paste0(out_dir,"plot_means_s1_",s1_label,".pdf"), p_mean, width=14, height=10)
  ggsave(paste0(out_dir,"plot_means_s1_",s1_label,".eps"), p_mean, width=14, height=10, device="eps")

  # ----------------------------------------------------------
  # Accumulate into surface
  # ----------------------------------------------------------

  surface_summary   <- rbind(surface_summary,   design_summary)
  surface_avg_sds   <- rbind(surface_avg_sds,   avg_sds_df)
  surface_avg_means <- rbind(surface_avg_means, avg_means_df)

} # end s1 loop

stopCluster(cl)

# ============================================================
# Save full combined surface
# ============================================================

save.image(file = paste0(master_dir, "surface_full.RData"))

write.csv(surface_summary,
          paste0(master_dir, "surface_design_summary.csv"), row.names=FALSE)
write.csv(surface_avg_sds,
          paste0(master_dir, "surface_avg_sds.csv"),        row.names=FALSE)
write.csv(surface_avg_means,
          paste0(master_dir, "surface_avg_means.csv"),      row.names=FALSE)

cat("\nSaved: surface_full.RData\n")
cat("Saved: surface_design_summary.csv\n")
cat("Saved: surface_avg_sds.csv\n")
cat("Saved: surface_avg_means.csv\n")


# ============================================================
# Raw optimal design points across full surface
# ============================================================

opt_tp    <- surface_summary[which.min(surface_summary$C1_tp),    ]
opt_logtp <- surface_summary[which.min(surface_summary$C1_logtp), ]

cat("\n--- Optimal design (C1_tp) ---\n");    print(opt_tp)
cat("\n--- Optimal design (C1_logtp) ---\n"); print(opt_logtp)

# ============================================================
# Kernel smoother — applied after loop using surface_summary
# ============================================================

smooth_1d <- function(tau_vec, c_vec, n_fine = 500) {
  h        <- tau_vec[2] - tau_vec[1]
  tau_fine <- seq(min(tau_vec), max(tau_vec), length.out = n_fine)
  c_smooth <- sapply(tau_fine, function(t) {
    w <- dnorm((tau_vec - t) / h)
    sum(w * c_vec) / sum(w)
  })
  opt_idx <- which.min(c_smooth)
  list(tau_fine = tau_fine,
       c_smooth = c_smooth,
       tau_star = tau_fine[opt_idx],
       c_star   = c_smooth[opt_idx])
}

# ============================================================
# Smoothed one-variable optimal designs — one row per s1
# ============================================================

cat("\n============================================================\n")
cat("Smoothed one-variable optimal designs\n")
cat("============================================================\n")

x1_vals_unique <- sort(unique(surface_summary$x1))

smoothed_rows <- list()
smooth_plot_data <- list()   # collect for per-s1 smooth plots

for (x1_val in x1_vals_unique) {

  sub      <- surface_summary[surface_summary$x1 == x1_val, ]
  sub      <- sub[order(sub$tau), ]
  s1_val   <- unique(sub$s1)
  s1_lab   <- as.character(round(1 / s1_val))

  sm_C1 <- smooth_1d(sub$tau, sub$C1_tp)
  sm_C2 <- smooth_1d(sub$tau, sub$C1_logtp)

  cat(sprintf("\n  s1 = 1/%s  (x1 = %.3f)\n", s1_lab, x1_val))
  cat(sprintf("    C1(tp)    : tau* = %.4f   value = %.6f\n", sm_C1$tau_star, sm_C1$c_star))
  cat(sprintf("    C2(logtp) : tau* = %.4f   value = %.6f\n", sm_C2$tau_star, sm_C2$c_star))

  smoothed_rows[[length(smoothed_rows) + 1]] <- data.frame(
    s1              = s1_val,
    x1              = x1_val,
    s1_inv          = round(1 / s1_val),
    tau_star_C1     = sm_C1$tau_star,
    C1_tp_smooth    = sm_C1$c_star,
    tau_star_C2     = sm_C2$tau_star,
    C1_logtp_smooth = sm_C2$c_star
  )

  # collect smooth plot data
  smooth_plot_data[[length(smooth_plot_data) + 1]] <- list(
    s1_lab  = s1_lab,
    x1_val  = x1_val,
    raw     = sub,
    sm_C1   = sm_C1,
    sm_C2   = sm_C2
  )
}

smoothed_onevariable <- do.call(rbind, smoothed_rows)

cat("\n--- Smoothed one-variable optimal designs (summary table) ---\n")
print(smoothed_onevariable)

write.csv(smoothed_onevariable,
          paste0(master_dir, "smoothed_onevariable_optimal.csv"),
          row.names = FALSE)
cat("Saved: smoothed_onevariable_optimal.csv\n")


# ============================================================
# Smoothed two-variable optimal design
# ============================================================

opt2_C1 <- smoothed_onevariable[which.min(smoothed_onevariable$C1_tp_smooth), ]
opt2_C2 <- smoothed_onevariable[which.min(smoothed_onevariable$C1_logtp_smooth), ]

cat("\n--- Smoothed two-variable optimal design (C1_tp) ---\n")
print(opt2_C1)

cat("\n--- Smoothed two-variable optimal design (C1_logtp) ---\n")
print(opt2_C2)

# ============================================================
# Per-s1 smooth plots: raw points + smooth curve
# saved to each s1 subfolder
# ============================================================

for (item in smooth_plot_data) {

  out_dir_s1 <- paste0(master_dir, "s1_", item$s1_lab, "/")

  smooth_df <- data.frame(
    tau       = c(item$sm_C1$tau_fine, item$sm_C2$tau_fine),
    value     = c(item$sm_C1$c_smooth, item$sm_C2$c_smooth),
    criterion = rep(c("C1(tp)", "C2(log tp)"), each = length(item$sm_C1$tau_fine))
  )

  raw_df <- data.frame(
    tau       = rep(item$raw$tau, 2),
    value     = c(item$raw$C1_tp, item$raw$C1_logtp),
    criterion = rep(c("C1(tp)", "C2(log tp)"), each = nrow(item$raw))
  )

  opt_df <- data.frame(
    criterion = c("C1(tp)", "C2(log tp)"),
    tau_star  = c(item$sm_C1$tau_star, item$sm_C2$tau_star),
    c_star    = c(item$sm_C1$c_star,   item$sm_C2$c_star)
  )

  p_smooth <- ggplot() +
    geom_point(data = raw_df,
               aes(x = tau, y = value),
               colour = "grey60", size = 1.5) +
    geom_line(data = smooth_df,
              aes(x = tau, y = value),
              colour = "black", linewidth = 0.9) +
    geom_vline(data = opt_df,
               aes(xintercept = tau_star),
               linetype = "dashed", colour = "red", linewidth = 0.6) +
    geom_point(data = opt_df,
               aes(x = tau_star, y = c_star),
               colour = "red", size = 3) +
    facet_wrap(~ criterion, scales = "free_y") +
    labs(x     = expression(tau),
         y     = "Criterion value",
         title = bquote("Smoothed criteria vs "*tau*
                        " | s1 = 1/"*.(item$s1_lab))) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(paste0(out_dir_s1, "plot_smooth_criteria_s1_", item$s1_lab, ".pdf"),
         p_smooth, width = 10, height = 5)
  ggsave(paste0(out_dir_s1, "plot_smooth_criteria_s1_", item$s1_lab, ".eps"),
         p_smooth, width = 10, height = 5, device = "eps")
  
  cat(sprintf("  Saved smooth plot for s1 = 1/%s\n", item$s1_lab))
}


# ============================================================
# 3D Surface plots — wireframe with 2D kernel-smoothed surface
# ============================================================

tau_vals <- sort(unique(surface_summary$tau))
x1_vals  <- sort(unique(surface_summary$x1))

h_tau <- tau_vals[2] - tau_vals[1]
h_x1  <- x1_vals[2]  - x1_vals[1]

tau_fine <- seq(min(tau_vals), max(tau_vals), length.out = 100)
x1_fine  <- seq(min(x1_vals),  max(x1_vals),  length.out = 50)

for (crit_col in c("C1_tp", "C1_logtp")) {
    
    z_smooth <- matrix(NA, nrow = length(tau_fine), ncol = length(x1_fine))
    
    for (ii in seq_along(tau_fine)) {
        for (jj in seq_along(x1_fine)) {
            w <- dnorm((surface_summary$tau - tau_fine[ii]) / h_tau) *
                dnorm((surface_summary$x1  - x1_fine[jj])  / h_x1)
            z_smooth[ii, jj] <- sum(w * surface_summary[[crit_col]]) / sum(w)
        }
    }
    
    opt_idx <- which(z_smooth == min(z_smooth), arr.ind = TRUE)
    cat(sprintf("\nSmoothed optimum for %s: tau=%.3f, x1=%.3f, value=%.4f\n",
                crit_col, tau_fine[opt_idx[1]], x1_fine[opt_idx[2]],
                z_smooth[opt_idx[1], opt_idx[2]]))
    
    surf_df   <- expand.grid(tau = tau_fine, x1 = x1_fine)
    surf_df$z <- as.vector(z_smooth)
    
    crit_lab <- ifelse(crit_col == "C1_tp",
                       expression(C[1](t[p])),
                       expression(C[1](log(t[p]))))
    
    # Crop to region of interest — shows convexity clearly
    x1_max_plot <- 0.5   # adjust this to taste — try 0.4 or 0.5
    x1_min_plot <- 0.1
    
    surf_df_crop <- subset(surf_df, x1 <= x1_max_plot & x1 >= x1_min_plot)
    
    p1 <- wireframe(z ~ tau * x1, data = surf_df_crop,
                    shade = FALSE, col = "grey70",
                    scales = list(arrows = FALSE),
                    xlab = expression(tau), ylab = expression(x[1]), zlab = "",
                    main = crit_lab, screen = list(z = 40, x = -65))
    
    p2 <- wireframe(z ~ tau * x1, data = surf_df_crop,
                    shade = FALSE, col = "grey70",
                    scales = list(arrows = FALSE),
                    xlab = expression(tau), ylab = expression(x[1]), zlab = "",
                    main = crit_lab, screen = list(z = 130, x = -65))
                    
    
    pdf(paste0(master_dir, "surface_3D_", crit_col, ".pdf"),
      width = 12, height = 6)
  print(p1, position = c(0, 0, 0.52, 1), more = TRUE)
  print(p2, position = c(0.48, 0, 1, 1), more = FALSE)
  dev.off()

  postscript(paste0(master_dir, "surface_3D_", crit_col, ".eps"),
             width = 12, height = 6,
             paper = "special", horizontal = FALSE)
  print(p1, position = c(0, 0, 0.52, 1), more = TRUE)
  print(p2, position = c(0.48, 0, 1, 1), more = FALSE)
  dev.off()

  cat(sprintf("Saved 3D surface: surface_3D_%s.pdf/.eps\n", crit_col))
}



# Reset par
par(mfrow = c(1,1))

# ============================================================
# Criterion curves: one curve per x1, plotted against tau
# ============================================================

for (crit_col in c("C1_tp","C1_logtp")) {

  p_curves <- ggplot(surface_summary,
                     aes(x     = tau,
                         y     = .data[[crit_col]],
                         colour = factor(round(x1,3)),
                         group  = factor(round(x1,3)))) +
    geom_line(linewidth=1.1) + geom_point(size=2.5) +
    scale_colour_viridis_d(name=expression(x[1])) +
    labs(x=expression(tau),
         y=crit_col,
         title=paste0(crit_col," vs tau — one curve per x1")) +
    theme_bw(base_size=14) + theme(plot.title=element_text(hjust=0.5))

  ggsave(paste0(master_dir,"surface_curves_",crit_col,".pdf"),
         p_curves, width=10, height=7)
  ggsave(paste0(master_dir,"surface_curves_",crit_col,".eps"),
         p_curves, width=10, height=7, device="eps")
  print(p_curves)
}

# ============================================================
# SD surface heatmaps: one per selected parameter
# (main parameters only — excludes theta to avoid noise)
# ============================================================

params_for_surface <- c("a1","a2","b1","b2","beta1","beta2","tp","log_tp")

for (param in params_for_surface) {

  col_name <- make.names(param)   # handles special chars like [
  if (!(col_name %in% names(surface_avg_sds))) next

  p_surf <- ggplot(surface_avg_sds,
                   aes(x=tau, y=x1, fill=.data[[col_name]])) +
    geom_tile() +
    scale_fill_viridis_c(name="Posterior SD", option="plasma") +
    labs(x=expression(tau), y=expression(x[1]),
         title=paste0("Posterior SD of ", param," — design surface")) +
    theme_bw(base_size=14) + theme(plot.title=element_text(hjust=0.5))

  fname <- gsub("[^a-zA-Z0-9]","_", param)
  ggsave(paste0(master_dir,"surface_SD_",fname,".pdf"),
         p_surf, width=9, height=6)
  ggsave(paste0(master_dir,"surface_SD_",fname,".eps"),
         p_surf, width=9, height=6, device="eps")
}

cat("\nAll outputs saved to:", master_dir, "\n")
cat("Done.\n")
