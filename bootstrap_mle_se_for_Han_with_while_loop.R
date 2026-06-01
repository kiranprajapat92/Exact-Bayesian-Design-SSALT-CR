rm(list = ls())  
library(parallel)
library(doParallel)
library(foreach)

set.seed(2026)

##########  Load functions and data #########

source("C:/Users/Kiran/Dropbox/RWTH/Newcastle_OBD_SSALT/all_functions.R")

solar_data <- data.frame(
  t = c(
    0.140, 0.783, 1.324, 1.582,
    1.716, 1.794, 1.883, 2.293,
    2.660, 2.674, 2.725, 3.085,
    3.924, 4.396, 4.612, 4.892,
    5.002, 5.022, 5.082, 5.112,
    5.147, 5.238, 5.244, 5.247,
    5.305, 5.337, 5.407, 5.408,
    5.445, 5.483, 5.717
  ),
  C_i = c(
    1,2,2,1, 2,2,2,2,
    2,2,2,2, 2,2,1,2,
    1,2,2,1, 1,1,1,1,
    1,2,1,2, 1,1,2
  )
)


#######  parameters #########
tc  <- 6
tau <- 5

s0 <- 1/293
s1 <- 1/293 #320.2136
s2 <- 1/353

x1 <- (s1 - s0) / (s2 - s0)
x2 <- (s2 - s0) / (s2 - s0)


##### Prepare likelihood inputs #######
solar_data <- solar_data[order(solar_data$t), ]

t_obs <- solar_data$t
C_obs <- solar_data$C_i

n  <- 35                       # planned sample size
n1 <- sum(t_obs <= tau)
nc <- sum(t_obs <= tc)

t_for_lik <- c(t_obs[t_obs <= tc], rep(tc, n - nc))
C_for_lik <- c(C_obs[t_obs <= tc], rep(NA, n - nc))

##### Original MLE #######
init_para <- c(log(3), -0.5, 1, log(3), -0.5, 1)

fit <- optim(
  par = init_para,
  fn  = likelihood_para_vec,
  t   = t_for_lik,
  C_i = C_for_lik,
  tau = tau,
  n1  = n1,
  nc  = nc,
  x1  = x1,
  x2  = x2,
  method = "L-BFGS-B",
  lower = c(-100,-100,1e-6,-100,-100,1e-6),
  upper = c(100,100,100,100,100,100)
)

para_hat <- fit$par
names(para_hat) <- c("a1","b1","beta1","a2","b2","beta2")
para_hat

# para_hat[2] <- -2

# t_q's
q=0.001
c(exp(para_hat[1]) * (-log(1 - q))^(1 / para_hat[3]), exp(para_hat[4]) * (-log(1 - q))^(1 / para_hat[6]));
#theta_ij's
c(exp(para_hat[1] + para_hat[2] * x1), exp(para_hat[4] + para_hat[5] * x1), exp(para_hat[1] + para_hat[2] * x2), exp(para_hat[4] + para_hat[5] * x2));

####  Parametric bootstrap for SEs ####
## Initialize

B_target <- 1000     # number of VALID bootstrap samples desired
b_valid  <- 0        # counter for valid samples
b_total  <- 0        # total attempts (for reporting)

boot_par <- matrix(NA, B_target, 6)
colnames(boot_par) <- c("a1","b1","beta1","a2","b2","beta2")

tq1_boot <- numeric(B_target)
tq2_boot <- numeric(B_target)

theta11_boot <- numeric(B_target)
theta12_boot <- numeric(B_target)
theta21_boot <- numeric(B_target)
theta22_boot <- numeric(B_target)

## bootstrap while loop

q <- 0.001
p_vec <- c(0.01, 0.1, 0.5)
tp_boot <- matrix(NA, B_target, length(p_vec))
colnames(tp_boot) <- paste0("p_", p_vec)

# set.seed(2026)

while (b_valid < B_target) {
  
  b_total <- b_total + 1
  
  ## 1. Simulate data
  sim <- simSSAT_CR(
    n = n,
    tau = tau,
    tc = tc,
    para_vec = para_hat,
    x1 = x1,
    x2 = x2
  )
  
  t_b <- pmin(sim$t, tc)
  C_b <- sim$C_i
  
  n1_b <- sum(t_b <= tau)
  nc_b <- sum(t_b <= tc)
  
  n_11 <- sum(t_b <= tau & C_b == 1)
  n_12 <- sum(t_b <= tau & C_b == 2)
  n_21 <- sum(t_b > tau & t_b <= tc & C_b == 1)
  n_22 <- sum(t_b > tau & t_b <= tc & C_b == 2)
  
  ## MLE existence condition
  if (n_11 == 0 || n_12 == 0 || n_21 == 0 || n_22 == 0)
    next
  
  ## 2. Refit model
  fit_b <- try(
    optim(
      par = para_hat,
      fn  = likelihood_para_vec,
      t   = t_b,
      C_i = C_b,
      tau = tau,
      n1  = n1_b,
      nc  = nc_b,
      x1  = x1,
      x2  = x2,
      method = "L-BFGS-B",
      lower = c(-100,-100,1e-6,-100,-100,1e-6),
      upper = c(100,100,100,100,100,100)
    ),
    silent = TRUE
  )
  
  if (inherits(fit_b, "try-error") || fit_b$convergence != 0)
    next
  
  ## 3. Store VALID replicate
  b_valid <- b_valid + 1
  
  boot_par[b_valid, ] <- fit_b$par
  
  tq1_boot[b_valid] <- exp(fit_b$par[1]) *
    (-log(1 - q))^(1 / fit_b$par[3])
  
  tq2_boot[b_valid] <- exp(fit_b$par[4]) *
    (-log(1 - q))^(1 / fit_b$par[6])
  
  theta11_boot[b_valid] <- exp(fit_b$par[1] + fit_b$par[2] * x1)
  theta12_boot[b_valid] <- exp(fit_b$par[4] + fit_b$par[5] * x1)
  theta21_boot[b_valid] <- exp(fit_b$par[1] + fit_b$par[2] * x2)
  theta22_boot[b_valid] <- exp(fit_b$par[4] + fit_b$par[5] * x2)
  
  a1_b <- fit_b$par[1]
  beta1_b <- fit_b$par[3]
  a2_b <- fit_b$par[4]
  beta2_b <- fit_b$par[6]
  
  theta1_b <- exp(a1_b)
  theta2_b <- exp(a2_b)
  
  for (k in 1:length(p_vec)) {
    
    p_now <- p_vec[k]
    
    f_root <- function(s) {
      log(1 - p_now) +
        (s/theta1_b)^beta1_b +
        (s/theta2_b)^beta2_b
    }
    
    tp_boot[b_valid, k] <- uniroot(f_root, lower = 1e-6, upper = 1000)$root
  }
  
}


cat("Valid bootstrap replicates :", B_target, "\n")
cat("Total bootstrap attempts   :", b_total, "\n")
cat("Acceptance rate            :",
    round(100 * B_target / b_total, 1), "%\n")

# boot_valid <- boot_par[complete.cases(boot_par), ]
# 
# ## Remove failed replications
# valid <- complete.cases(boot_par)
# 
# boot_par   <- boot_par[valid, ]
# tq1_boot   <- tq1_boot[valid]
# tq2_boot   <- tq2_boot[valid]
# 
# theta11_boot <- theta11_boot[valid]
# theta12_boot <- theta12_boot[valid]
# theta21_boot <- theta21_boot[valid]
# theta22_boot <- theta22_boot[valid]

## Bootstrap means and SEs
par_mean <- apply(boot_par, 2, mean)
par_se   <- apply(boot_par, 2, sd)

tp_mean <- apply(tp_boot, 2, mean)
tp_se   <- apply(tp_boot, 2, sd)

pth_quantile_summary <- data.frame(
  p    = p_vec,
  Mean = tp_mean,
  SE   = tp_se
)

pth_quantile_summary


quantile_table <- data.frame(
  Cause = c("1", "2"),
  q     = q,
  Mean  = c(mean(tq1_boot), mean(tq2_boot)),
  SE    = c(sd(tq1_boot),   sd(tq2_boot))
)

theta_table <- data.frame(
  Parameter = c("theta11","theta12","theta21","theta22"),
  Mean = c(
    mean(theta11_boot),
    mean(theta12_boot),
    mean(theta21_boot),
    mean(theta22_boot)
  ),
  SE = c(
    sd(theta11_boot),
    sd(theta12_boot),
    sd(theta21_boot),
    sd(theta22_boot)
  )
)

par_mean       # para_vaec
par_se
quantile_table # quantile for j-the cause
theta_table    # scale parameters


## Bootstrap CIs

boot_ci_percentile <- function(x, alpha = 0.05) {
  quantile(x, probs = c(alpha/2, 1 - alpha/2))
}

## para_vec
par_ci_perc <- t(apply(
  boot_par, 2, boot_ci_percentile
))

colnames(par_ci_perc) <- c("Lower", "Upper")
par_ci_perc


tp_ci <- t(apply(tp_boot, 2, boot_ci_percentile))
colnames(tp_ci) <- c("Lower", "Upper")

quantile_ci_table <- cbind(
  p = p_vec,
  tp_ci
)

quantile_ci_table


##  Weibull quantiles for j-th cause
tq1_ci_perc <- boot_ci_percentile(tq1_boot)
tq2_ci_perc <- boot_ci_percentile(tq2_boot)

tq_ci_table <- data.frame(
  Cause = c("1","2"),
  q = q,
  Lower = c(tq1_ci_perc[1], tq2_ci_perc[1]),
  Upper = c(tq1_ci_perc[2], tq2_ci_perc[2])
)

tq_ci_table

## scale parameters

theta_ci_perc <- rbind(
  theta11 = boot_ci_percentile(theta11_boot),
  theta12 = boot_ci_percentile(theta12_boot),
  theta21 = boot_ci_percentile(theta21_boot),
  theta22 = boot_ci_percentile(theta22_boot)
)

colnames(theta_ci_perc) <- c("Lower", "Upper")
theta_ci_perc

# Better visible results:
par_results <- cbind(par_mean,par_se,par_ci_perc); par_results
pth_quantile_results <- cbind(pth_quantile_summary,quantile_ci_table)[,-4];pth_quantile_results;
quantile_results <- cbind(quantile_table, tq_ci_table)[,-c(5,6)]; quantile_results
theta_results <- cbind(theta_table, theta_ci_perc)[,-1]; theta_results


# par(mfrow = c(2,3), mar = c(4,4,1,1))
# hist(boot_par[,1], probability=TRUE, breaks=25,
#      xlab = expression(a[1]), main = "")
# lines(density(boot_par[,1]), lwd=2)
# hist(boot_par[,2], probability=TRUE, breaks=25,
#      xlab = expression(b[1]), main = "")
# lines(density(boot_par[,2]), lwd=2)
# hist(boot_par[,3], probability=TRUE, breaks=30,
#      xlab = expression(beta[1]), main = "")
# lines(density(boot_par[,3]), lwd=2)
# hist(boot_par[,4], probability=TRUE, breaks=25,
#      xlab = expression(a[2]), main = "")
# lines(density(boot_par[,4]), lwd=2)
# hist(boot_par[,5], probability=TRUE, breaks=20,
#      xlab = expression(b[2]), main = "")
# lines(density(boot_par[,5]), lwd=2)
# hist(boot_par[,6], probability=TRUE, breaks=20,
#      xlab = expression(beta[2]), main = "")
# lines(density(boot_par[,6]), lwd=2)
# 
# 
# par(mfrow = c(2,3), mar = c(4,4,1,1))
# hist(theta11_boot, probability=TRUE, breaks=40,
#      xlab = expression(theta[11]), main = "")
# lines(density(theta11_boot), lwd=2)
# hist(theta12_boot, probability=TRUE, breaks=45,
#      xlab = expression(theta[12]), main = "")
# lines(density(theta12_boot), lwd=2)
# hist(theta21_boot, probability=TRUE, breaks=40,
#      xlab = expression(theta[21]), main = "")
# lines(density(theta21_boot), lwd=2)
# hist(theta22_boot, probability=TRUE, breaks=40,
#      xlab = expression(theta[22]), main = "")
# lines(density(theta22_boot), lwd=2)
# hist(tq1_boot, probability=TRUE, breaks=45,
#      xlab = expression(t["0,1"]^q), main = "")
# lines(density(tq1_boot), lwd=2)
# hist(tq2_boot, probability=TRUE, breaks=45,
#      xlab = expression(t["0,2"]^q), main = "")
# lines(density(tq2_boot), lwd=2)
# 
# par(mfrow = c(1,3), mar = c(4,4,1,1))
# hist(tp_boot[,1], probability=TRUE, breaks=20,
#      xlab = expression(t[0.01](x["0"])), main = "")
# lines(density(tp_boot[,1]), lwd=2) 
# hist(tp_boot[,2], probability=TRUE, breaks=15,
#      xlab = expression(t[0.1](x["0"])), main = "")
# lines(density(tp_boot[,2]), lwd=2)
# hist(tp_boot[,3], probability=TRUE, breaks=20,
#      xlab = expression(t[0.5](x["0"])), main = "")
# lines(density(tp_boot[,3]), lwd=2)
# 
# 

# ============================================================
# 3x3 Bootstrap histogram plots
# ============================================================

library(ggplot2)
library(gridExtra)

hist_colour <- "#1B2A7A"
dens_colour <- "#8B0000"
# hist_colour <- "#4A7FB5"
# dens_colour <- "#1A1A2E"

plot_hist <- function(x, xlab, bw_adjust = 0.6) {
  df <- data.frame(x = x)
  ggplot(df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)),
                   fill    = hist_colour,
                   colour  = hist_colour,
                   alpha   = 0.75,
                   bins    = 40) +
    geom_density(colour    = dens_colour,
                 linewidth = 0.9,
                 adjust    = bw_adjust) +
    labs(x = xlab, y = NULL) +
    theme_bw(base_size = 13) #+
    #theme(axis.text.y      = element_blank(),
    #      axis.ticks.y     = element_blank(),
    #      panel.grid.minor = element_blank())
}

# ── Row 1: cause 1 reparametrised parameters ─────────────
p1 <- plot_hist(tq1_boot,
                expression(varphi[11] == t[paste("0,1")]^q),
                bw_adjust = 0.3)
p2 <- plot_hist(-boot_par[, 2],
                expression(varphi[21] == -b[1]),
                bw_adjust = 0.8)
p3 <- plot_hist(boot_par[, 3],
                expression(varphi[31] == beta[1]),
                bw_adjust = 0.3)

# ── Row 2: cause 2 reparametrised parameters ─────────────
p4 <- plot_hist(tq2_boot,
                expression(varphi[12] == t[paste("0,2")]^q),
                bw_adjust = 0.5)
p5 <- plot_hist(-boot_par[, 5],
                expression(varphi[22] == -b[2]),
                bw_adjust = 0.8)
p6 <- plot_hist(boot_par[, 6],
                expression(varphi[32] == beta[2]),
                bw_adjust = 0.8)

# ── Row 3: quantiles at use stress ───────────────────────
p7 <- plot_hist(tp_boot[, 1],
                expression(t[0.01](x[0])),
                bw_adjust = 0.6)
p8 <- plot_hist(tp_boot[, 2],
                expression(t[0.10](x[0])),
                bw_adjust = 0.8)
p9 <- plot_hist(tp_boot[, 3],
                expression(t[0.50](x[0])),
                bw_adjust = 0.8)

# ── Arrange and save ──────────────────────────────────────
fig_all <- grid.arrange(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3
)

ggsave("C:/Users/Kiran/Desktop/plots/bootstrap_hist_all.pdf",
       fig_all, width = 12, height = 10)
ggsave("C:/Users/Kiran/Desktop/plots/bootstrap_hist_all.eps",
       fig_all, width = 12, height = 10, device = "eps")