
# ---------- simulate observations under simple SSALT with competeing risks ----------

simSSAT_CR <- function(n, tau, tc, para_vec, x1, x2)
{
  a1 <- para_vec[1]; b1 <- para_vec[2]
  a2 <- para_vec[4]; b2 <- para_vec[5]
  betaa <- c(para_vec[3], para_vec[6])
  thetaa <- matrix(0, nrow = 2,ncol = 2)
  thetaa[1,1] <- exp(a1 + b1 * x1)
  thetaa[1,2] <- exp(a2 + b2 * x1)
  thetaa[2,1] <- exp(a1 + b1 * x2)
  thetaa[2,2] <- exp(a2 + b2 * x2)
  
  z <- rep(0, n)
  for(j in 1:2)
  {
    zj <- rweibull(n, betaa[j], thetaa[1,j])
    z <- cbind(z, zj)
  }
  z <- z[,-1]
  t1 <- apply(z,1,min)
  get_min_column_index <- function(row) { which.min(row) }
  cause_of_failure <- apply(z, 1, get_min_column_index)
  sorted_indices <- order(t1)
  t1 <- t1[sorted_indices]
  cause_of_failure <- cause_of_failure[sorted_indices]
  data1 <- cbind(t1, cause_of_failure)
  z1 <- t1[t1 <= tau]
  n1 <- length(z1)
  data1  <- if(n1 >= 2){
    data1[1:n1,]
  }else{
    if(n1==1){
      data1[1:n1,]
    }else{
      numeric()
    }
  }
  n2 <- n-n1
  if(n2>0){
    z2 <- rep(0, n2)
    for(j in 1:2)
    {
      uj <- runif(n2)
      z2j<- thetaa[2,j] * (-log(uj*exp(-(tau/thetaa[1,j])^betaa[j])))^(1/betaa[j]) - tau * ((thetaa[2,j]/thetaa[1,j])-1)
      z2 <- cbind(z2, z2j)
    }
    z2 <- z2[,-1]
    t2 <- if(n2==1){
      min(z2)
    }else{
      apply(z2,1,min)
    }
    cause_of_failure <- if(n2==1){
      which.min(array(z2)) 
    }else{
      apply(z2, 1, get_min_column_index)
    }
    
    sorted_indices <- order(t2)
    t2 <- t2[sorted_indices]
    cause_of_failure <- cause_of_failure[sorted_indices]
    data2 <- cbind(t2, cause_of_failure) 
  }
  data <- if(n2>0){rbind(data1, data2)}else{data1}
  list(t = data[,1], C_i = data[,2])
}



# ---------- Likelihood function ------------


likelihood_para_vec <- function(para_vec, t, C_i, tau, n1, nc, x1, x2) {
  a1 <- para_vec[1]; b1 <- para_vec[2]; beta1 <- para_vec[3]
  a2 <- para_vec[4]; b2 <- para_vec[5]; beta2 <- para_vec[6]
  thetaa <- matrix(NA, nrow = 2, ncol = 2)
  thetaa[1,1] <- exp(a1 + b1 * x1)
  thetaa[1,2] <- exp(a2 + b2 * x1)
  thetaa[2,1] <- exp(a1 + b1 * x2)
  thetaa[2,2] <- exp(a2 + b2 * x2)
  sai_l1 <- function(j1, t, tau, thetaa) { (1 / thetaa[1, j1]) * t }
  sai_l2 <- function(j1, t, tau, thetaa) { ((t - tau)/thetaa[2, j1]) + (tau/thetaa[1, j1]) }
  n <- length(t)
  log_L <- rep(NA, n); 
  # log_lik <- numeric(n)
  sai_l11 <- sai_l1(1, t, tau, thetaa); sai_l12 <- sai_l1(2, t, tau, thetaa)
  sai_l21 <- sai_l2(1, t, tau, thetaa); sai_l22 <- sai_l2(2, t, tau, thetaa)
  if (n1 >= 1) {
    for (i in 1:n1) {
      if (C_i[i] == 1) {
        log_L[i] <- log(beta1) - log(thetaa[1,1]) + (beta1 - 1) * log(sai_l11[i]) -
          sai_l11[i]^beta1 - sai_l12[i]^beta2
      } else {
        log_L[i] <- log(beta2) - log(thetaa[1,2]) + (beta2 - 1) * log(sai_l12[i]) -
          sai_l11[i]^beta1 - sai_l12[i]^beta2
      }
      # log_lik[i] <- exp(log_L[i])
      # if (log_lik[i] < 1e-10) log_lik[i] <- 1e-10
    }
  }
  if (nc >= (n1+1)) {
    for (i in (n1+1):nc) {
      if (C_i[i] == 1) {
        log_L[i] <- log(beta1) - log(thetaa[2,1]) + (beta1 - 1) * log(sai_l21[i]) -
          sai_l21[i]^beta1 - sai_l22[i]^beta2
      } else {
        log_L[i] <- log(beta2) - log(thetaa[2,2]) + (beta2 - 1) * log(sai_l22[i]) -
          sai_l21[i]^beta1 - sai_l22[i]^beta2
      }
      # log_lik[i] <- exp(log_L[i])
      # if (log_lik[i] < 1e-9) log_lik[i] <- 1e-9
    }
  }
  if (nc < n) {
    for (i in (nc+1):n) {
      log_L[i] <- - (sai_l21[i]^beta1) - (sai_l22[i]^beta2)
      # log_lik[i] <- exp(log_L[i])
      # if (log_lik[i] < 1e-9) log_lik[i] <- 1e-9
    }
  }
  # prob <- sum(log(log_lik))
  prob <- sum(log_L)
  return(-prob)
}



# ---------- CDF (not censored or truncated) ------------


cdf_stepstress <- function(t, tau, thetaa, beta1, beta2) {
  
  Lambda <- numeric(length(t))
  
  for (i in seq_along(t)) {
    if (t[i] <= tau) {
      Lambda[i] <- (t[i]/thetaa[1,1])^beta1 +
        (t[i]/thetaa[1,2])^beta2
    } else {
      Lambda[i] <- (tau/thetaa[1,1] + (t[i]-tau)/thetaa[2,1])^beta1 +
        (tau/thetaa[1,2] + (t[i]-tau)/thetaa[2,2])^beta2
    }
  }
  
  return(1 - exp(-Lambda))
}


# ----------  Empirical CDF --------------


empirical_cdf <- function(t, tc) {
  t <- sort(t)
  n <- length(t_for_lik)
  F_i <- 0
  for (i in 1:length(t)) {
    F_i <- cbind(F_i,sum(t_for_lik<=t[i])/n)
  }
  return(F_i[-1])
}


# ----------  statistics value --------------


CvM_stat_fun <- function(t_sample) {
  t_sample <- sort(t_sample)
  F_modelCvM <- cdf_stepstress(
    t_sample, tau, thetaa_hat,
    para_hat[3], para_hat[6]
  )
  
  BB <- 1 / (12*n) 
  for(ind3 in 1:n)
  {
    BB <- BB + (((2*ind3-1)/(2*n)) - F_modelCvM[ind3])^2
  }
  BB
}
