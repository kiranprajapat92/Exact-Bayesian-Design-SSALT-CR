rm(list = ls())         # to remove the memory related to previously assigned variables
library(parallel)
library(doParallel)
library(cmdstanr)
library(bayesplot)
library(ggplot2)
par(cex = 1.25, mar = c(4.75, 4.85, 2.25, 1.5))

source("C:/Users/Kiran/Dropbox/RWTH/Newcastle_OBD_SSALT/all_functions.R")

strt <- Sys.time()

set.seed(2026) #135

###########################################################################

# ============================================================
# PRIOR SELECTION (CHANGE ONLY THIS)
# ============================================================
# 1 = Uniform (mean ± 2SE)
# 2 = Uniform (mean ± 3SE)
# 3 = Uniform (bootstrap percentile CI)
# 4 = Gamma (mean/SE matched)
# 5 = Normal (mean/SE matched)


# prior_type <- 4   # 1,2,3,4,5

#######  parameters #########

n = 35
tc  <- 6
tau <- 5

s0 <- 1/293
s1 <- 1/320.2136
s2 <- 1/353

x0 <- 0
x1 <- (s1 - s0) / (s2 - s0)
x2 <- (s2 - s0) / (s2 - s0)

para_hat <- c(4.5064079, -4.7131110,  0.7692292,  2.0409831, -1.2277314,  1.5320872)

para_true <- para_hat



## Gamma 
# Cause 1
hpar11 <- c(0.1634133, 0.3704860) # t_{0,1}^q
hpar21 <- c(4.280502, 1.2736822)  # -b1
hpar31 <- c(1.200586, 1.2724302)  # beta1

# Cause 2
hpar12 <- c(0.1526861, 0.1549505) # t_{0,2}^q
hpar22 <- c(1.402464, 0.5039255)  # -b2
hpar32 <- c(1.698900, 0.4604109)  # beta2

hpars <- list(
  cause1 = list(hpar11, hpar21, hpar31),
  cause2 = list(hpar12, hpar22, hpar32)
)

gamma_hyper <- function(mean_sd) {
  mean_x <- mean_sd[1]
  sd_x   <- mean_sd[2]
  shape  <- mean_x^2 / sd_x^2
  rate   <- mean_x / sd_x^2
  c(shape, rate)
}

gphi11 <- gamma_hyper(hpar11)
gphi21 <- gamma_hyper(hpar21)
gphi31 <- gamma_hyper(hpar31)

gphi12 <- gamma_hyper(hpar12)
gphi22 <- gamma_hyper(hpar22)
gphi32 <- gamma_hyper(hpar32)


##############################################################

a1 <- para_true[1]
b1 <- para_true[2]
betaa1 <- para_true[3]
a2 <- para_true[4]
b2 <- para_true[5]
betaa2 <- para_true[6]
# Quantile of lifetime distribution under SSALT
p <- 0.01
funcc_fs <- function(s)
{
  f_s <- log(1 - p) + (s * exp(-a1))^betaa1 + (s * exp(-a2))^betaa2
  f_s
}
yx <- 0
xx <- seq(0.0001,10,by = 0.1)
for(i in 1:length(xx))
{
  yx <- c(yx, funcc_fs(xx[i]))
}
yx <- yx[-1]
plot(xx,yx, ylab = "f_s", xlab="s", "l")
lines(xx,rep(0,length(xx)))

tp = 5 
for (ii in 1:1000) {
  f_s = log(1 - p) + (tp * exp(-a1))^betaa1 + (tp * exp(-a2))^betaa2
  print(ii)
  print(f_s)
  if (abs(f_s) < 0.0001) {
    break  # Stop the loop if f_s is close to zero
  }
  tp = tp - f_s / ((betaa1 * (tp^(betaa1-1)) * ((exp(-a1))^betaa1)) + (betaa2 * (tp^(betaa2-1)) * ((exp(-a2))^betaa2)))
}
tp_x0 <- tp;
tp_x0


# Transformed parameters for priors
q <- 0.001
tq_1 <- exp(para_hat[1]) * (-log(1 - q))^(1 / para_hat[3]);
tq_2 <- exp(para_hat[4]) * (-log(1 - q))^(1 / para_hat[6]);
#theta_ij's
matrix(c(exp(para_hat[1] + para_hat[2] * x1), exp(para_hat[4] + para_hat[5] * x1), exp(para_hat[1] + para_hat[2] * x2), exp(para_hat[4] + para_hat[5] * x2)), byrow = TRUE,nrow = 2);

trans_para_cause1 <- c(tq_1, -b1, betaa1)
trans_para_cause2 <- c(tq_1, -b2, betaa2)
round(as.data.frame(matrix(c(trans_para_cause1,trans_para_cause2),2,3,byrow=TRUE), row.names = c("cause1","cause2")),digits = 3)

###########################################################################

## 1. Simulate data
sim <- simSSAT_CR(
  n = n,
  tau = tau,
  tc = tc,
  para_vec = para_true,
  x1 = x1,
  x2 = x2
)

t <- pmin(sim$t, tc)
C_i <- sim$C_i

n1<-sum(t<=tau)
nc<-sum(t<=tc)
############################################################


file <- file.path(cmdstan_path(), "examples/bernoulli/bernoulli_han_ssalt_gamma_repara.stan")
options("cmdstanr_verbose"=TRUE)
mod <- cmdstan_model(file, compile = TRUE,
                     cpp_options = list(stan_threads = TRUE))


init_fun <- function() {
  list(
    phi_11 = max(rgamma(1, shape = gphi11[1], rate = gphi11[2]), 1e-2),
    phi_21 = max(rgamma(1, shape = gphi21[1], rate = gphi21[2]), 1e-2),
    phi_31 = max(rgamma(1, shape = gphi31[1], rate = gphi31[2]), 1e-2),
    phi_12 = max(rgamma(1, shape = gphi12[1], rate = gphi12[2]), 1e-2),
    phi_22 = max(rgamma(1, shape = gphi22[1], rate = gphi22[2]), 1e-2),
    phi_32 = max(rgamma(1, shape = gphi32[1], rate = gphi32[2]), 1e-2)
  )
}


parameters<-c("a1", "a2", "b1", "b2","beta1", "beta2", "tp", "log_tp")

data_list <- list(
  n = n,
  tau = tau,
  tc = tc,
  n1 = n1,
  nc = nc,
  t = t,
  C_i = C_i,
  x1 = x1,
  x2 = x2,
  x0 = x0,
  p = p,
  
  # lognormal priors
  gphi11 = gphi11,
  gphi21 = gphi21,
  gphi31 = gphi31,
  gphi12 = gphi12,
  gphi22 = gphi22,
  gphi32 = gphi32
)


json_file <- tempfile(fileext = ".json")
write_stan_json(data_list, json_file)
# cat(readLines(json_file), sep = "\n")

SSAT.sim <- mod$sample(
  data = json_file,
  chains = 3,
  parallel_chains = 3,
  init = init_fun,
  iter_warmup = 1000,
  threads_per_chain = 1,
  iter_sampling = 1000, 
  # thin = 5,
  # seed = 135,
  show_exceptions = TRUE
)
# Extracting summary statistics
# SSAT.sim$print(max_rows = 20)
tp_log_tp_summary <- SSAT.sim$summary(variables = c("tp", "log_tp"), "sd")
sd_tp <- tp_log_tp_summary[1,2]
sd_log_tp <- tp_log_tp_summary[2,2]

########## Fit Check ############

check_stan_fit <- function(fit) {
  
  summ <- fit$summary()
  
  # 1. Rhat check
  rhat_ok <- all(summ$rhat < 1.01, na.rm = TRUE)
  
  # 2. Divergences
  sampler_diag <- fit$sampler_diagnostics()
  n_div <- sum(sampler_diag[, , "divergent__"])
  div_ok <- (n_div == 0)
  
  # 3. Treedepth
  n_treedepth <- sum(sampler_diag[, , "treedepth__"] >= 
                       fit$metadata()$max_treedepth)
  tree_ok <- (n_treedepth == 0)
  
  list(
    rhat_ok = rhat_ok,
    div_ok = div_ok,
    tree_ok = tree_ok,
    all_ok = rhat_ok & div_ok & tree_ok,
    n_divergent = n_div,
    n_treedepth = n_treedepth
  )
}

check_stan_fit(SSAT.sim)
#############################################################


# color_scheme_set("darkgray")
# 
# my_theme <- ggplot2::theme_gray() + theme(
#   legend.text  = element_text(size = 16),
#   legend.title = element_text(size = 16),
#   axis.text    = element_text(size = 14),
#   axis.title   = element_text(size = 15),
#   strip.text   = element_text(size = 16)
# )
# 
# # Plot 1: a1, a2, b1, b2
# # Export -> Save as EPS -> Width: 800, Height: 700
# 
# my_labeller1 <- as_labeller(
#   x = c(
#     'a1' = 'a[1]',
#     'a2' = 'a[2]',
#     'b1' = 'b[1]',
#     'b2' = 'b[2]'
#   ),
#   default = label_parsed
# )
# 
# p1 <- mcmc_combo(
#   x       = SSAT.sim$draws(),
#   combo   = c("dens_overlay", "trace"),
#   pars    = parameters[1:4],
#   facet_args = list(
#     ncol     = 1,
#     labeller = my_labeller1
#   ),
#   gg_theme = my_theme
# )
# print(p1)
# 
# 
# # Plot 2: beta1, beta2, tp, log_tp
# # Export -> Save as EPS -> Width: 800, Height: 700
# 
# my_labeller2 <- as_labeller(
#   x = c(
#     'beta1'  = 'beta[1]',
#     'beta2'  = 'beta[2]',
#     'tp'     = 't[p](x[0])',
#     'log_tp' = 'log(t[p](x[0]))'
#   ),
#   default = label_parsed
# )
# 
# p2 <- mcmc_combo(
#   x       = SSAT.sim$draws(),
#   combo   = c("dens_overlay", "trace"),
#   pars    = parameters[5:8],
#   facet_args = list(
#     ncol     = NULL,
#     labeller = my_labeller2
#   ),
#   gg_theme = my_theme
# )
# print(p2)
# 
# 
# # Plot 3: ACF
# # Export -> Save as PDF -> Width: 12, Height: 8
# 
# my_labeller3 <- as_labeller(
#   x = c(
#     'a1'     = 'a[1]',
#     'a2'     = 'a[2]',
#     'b1'     = 'b[1]',
#     'b2'     = 'b[2]',
#     'beta1'  = 'beta[1]',
#     'beta2'  = 'beta[2]',
#     'tp'     = 't[p](x[0])',
#     'log_tp' = 'log(t[p](x[0]))',
#     '1'      = '1',
#     '2'      = '2',
#     '3'      = '3'
#   ),
#   default = label_parsed
# )
# 
# p3 <- mcmc_acf(
#   SSAT.sim$draws(),
#   pars = parameters,
#   lags = 15,
#   facet_args = list(
#     labeller = my_labeller3
#   )
# ) + theme(
#   strip.text = element_text(size = 17),
#   axis.text  = element_text(size = 15),
#   axis.title = element_text(size = 16)
# )
# print(p3)
# 
# # # mcmc_dens(SSAT.sim$draws(), pars = parameters)
# 
# mcmc_chains <- SSAT.sim$draws()
# 
# # Choose the relevant columns for output
# output_data <- SSAT.sim$summary()
# # Specify the output file path
# output_file <- "C:/Users/Kiran Prajapat/OneDrive/xiao_liu_data/diago_temp2/cmdstan_output.csv"
# # # Save the data.frame to a CSV file
# # write.csv(output_data, file = output_file, row.names = FALSE)
# 
# SSAT.sim$print(max_rows = 20)
# SSAT.sim$summary(variables = c("a1", "b1", "beta1", "a2", "b2", "beta2","tp", "log_tp"), c("mean","sd"))
# 
# # round(as.data.frame(matrix(c(para_cause1,para_cause2),2,5,byrow=TRUE), row.names = c("cause1","cause2")),digits = 3)
# # round(as.data.frame(matrix(c(trans_para_cause1,trans_para_cause2),2,3,byrow=TRUE), row.names = c("cause1","cause2")),digits = 3)

library(ggplot2)
library(gridExtra)

# ---- base theme ----
base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    plot.title       = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.title       = element_text(size = 11),
    axis.text        = element_text(size = 8),
    legend.position  = "none",
    plot.margin      = margin(4, 6, 4, 6)
  )

# ---- darker, richer base colours ----
base_colours <- c(
  a1     = "#1A4A5E",
  a2     = "#A8115A",
  b1     = "#C99D00",
  b2     = "#9B1C1C",
  beta1  = "#0D57A1",
  beta2  = "#4A148C",
  tp     = "#BF2000",
  log_tp = "#001064"
)

# ---- 3 shades: dark, medium, light ----
make_shades <- function(base_col) {
  light  <- colorRampPalette(c(base_col, "#FFFFFF"))(5)
  shades <- c(base_col, light[1], light[2])
  names(shades) <- c("1", "2", "3")
  shades
}


# ---- labels ----
lab <- list(
  a1     = "a[1]",
  a2     = "a[2]",
  b1     = "b[1]",
  b2     = "b[2]",
  beta1  = "beta[1]",
  beta2  = "beta[2]",
  tp     = "t[p](x[0])",
  log_tp = "log(t[p](x[0]))"
)

# ---- main plotting function ----
plot_trace_dens_acf <- function(fit, pars, labeller_list,
                                base_colours, n_lags = 15) {
  plot_list <- list()
  
  for (k in seq_along(pars)) {
    
    par      <- pars[k]
    base_col <- base_colours[par]
    shades   <- make_shades(base_col)
    par_lab  <- labeller_list[[par]]
    
    draws_arr <- fit$draws(variables = par)
    n_iter    <- dim(draws_arr)[1]
    n_chains  <- dim(draws_arr)[2]
    
    # long format — explicit chain extraction
    df <- do.call(rbind, lapply(1:n_chains, function(ch) {
      data.frame(
        iteration = 1:n_iter,
        value     = as.vector(draws_arr[, ch, 1]),
        chain     = factor(ch, levels = c("1","2","3"))
      )
    }))
    
    # --- Trace ---
    post_mean <- mean(as.vector(draws_arr[, , 1]))
    
    p_trace <- ggplot(df, aes(x = iteration, y = value,
                              colour = chain, group = chain)) +
      geom_line(linewidth = 0.25, alpha = 0.95) +
      scale_colour_manual(values = shades,
                          breaks = c("1","2","3")) +
      geom_hline(yintercept = post_mean,
                 linetype   = "dashed",
                 colour     = "black",
                 linewidth  = 0.7) +
      labs(title = "Trace Plot",
           x     = "Iteration",
           y     = parse(text = par_lab)) +
      base_theme
    
    # --- Density ---
    p_dens <- ggplot(df, aes(x = value,
                             colour = chain,
                             # fill   = chain,
                             group  = chain)) +
      geom_density(alpha = 0.30, linewidth = 0.7) +
      scale_colour_manual(values = shades,
                          breaks = c("1","2","3")) +
      scale_fill_manual(values   = shades,
                        breaks   = c("1","2","3")) +
      labs(title = "Posterior Density",
           x     = parse(text = par_lab),
           y     = "Density") +
      base_theme
    
    # --- ACF ---
    acf_list <- lapply(1:n_chains, function(ch) {
      vals     <- as.vector(draws_arr[, ch, 1])
      acf_vals <- acf(vals, lag.max = n_lags, plot = FALSE)$acf[-1]
      data.frame(
        lag   = 1:n_lags,
        acf   = acf_vals,
        chain = factor(ch, levels = c("1","2","3"))
      )
    })
    df_acf <- do.call(rbind, acf_list)
    ci     <- qnorm(0.975) / sqrt(n_iter)
    
    p_acf <- ggplot(df_acf, aes(x = lag, y = acf,
                                fill = chain, group = chain)) +
      geom_bar(stat     = "identity",
               position = position_dodge(width = 0.8),
               width    = 0.7,
               alpha    = 0.9) +
      scale_fill_manual(values = shades,
                        breaks = c("1","2","3")) +
      geom_hline(yintercept =  ci,
                 linetype = "dashed", colour = "blue",
                 linewidth = 0.5) +
      geom_hline(yintercept = -ci,
                 linetype = "dashed", colour = "blue",
                 linewidth = 0.5) +
      geom_hline(yintercept = 0,
                 colour = "black", linewidth = 0.3) +
      scale_x_continuous(breaks = seq(0, n_lags, by = 5)) +
      scale_y_continuous(limits = c(-0.25, 1.0),
                         breaks = seq(-0.2, 1.0, by = 0.2)) +
      labs(title = "ACF",
           x     = "Lag",
           y     = "Autocorrelation") +
      base_theme
    
    plot_list <- c(plot_list, list(p_trace, p_dens, p_acf))
  }
  
  gridExtra::grid.arrange(
    grobs = plot_list,
    ncol  = 3,
    nrow  = length(pars)
  )
}

# ---- run and save ----

# Figure 1: a1, a2, b1, b2
p_fig1 <- plot_trace_dens_acf(
  fit           = SSAT.sim,
  pars          = c("a1","a2","b1","b2"),
  labeller_list = lab,
  base_colours  = base_colours,
  n_lags        = 15
)

# Figure 2: beta1, beta2, tp, log_tp
p_fig2 <- plot_trace_dens_acf(
  fit           = SSAT.sim,
  pars          = c("beta1","beta2","tp","log_tp"),
  labeller_list = lab,
  base_colours  = base_colours,
  n_lags        = 15
)

# save as EPS
postscript("fig_diag_ab.eps",
           width = 15, height = 14,
           paper = "special", horizontal = FALSE)
plot_trace_dens_acf(SSAT.sim,
                    pars          = c("a1","a2","b1","b2"),
                    labeller_list = lab,
                    base_colours  = base_colours,
                    n_lags        = 15)
dev.off()

postscript("fig_diag_betatp.eps",
           width = 15, height = 14,
           paper = "special", horizontal = FALSE)
plot_trace_dens_acf(SSAT.sim,
                    pars          = c("beta1","beta2","tp","log_tp"),
                    labeller_list = lab,
                    base_colours  = base_colours,
                    n_lags        = 15)
dev.off()

# also save as PDF
pdf("fig_diag_ab.pdf", width = 15, height = 14)
plot_trace_dens_acf(SSAT.sim,
                    pars          = c("a1","a2","b1","b2"),
                    labeller_list = lab,
                    base_colours  = base_colours,
                    n_lags        = 15)
dev.off()

pdf("fig_diag_betatp.pdf", width = 15, height = 14)
plot_trace_dens_acf(SSAT.sim,
                    pars          = c("beta1","beta2","tp","log_tp"),
                    labeller_list = lab,
                    base_colours  = base_colours,
                    n_lags        = 15)
dev.off()
