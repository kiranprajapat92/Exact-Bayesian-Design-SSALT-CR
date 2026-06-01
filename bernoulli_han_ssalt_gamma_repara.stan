functions {

  // -------------------------------------------------------
  // Reparameterised: accept log_thetaa instead of thetaa
  // so all internal divisions use exp(log_theta) which
  // is always finite and well-behaved.
  // -------------------------------------------------------

  vector sai_l1(int j, vector t, matrix log_thetaa) {
    return t * exp(-log_thetaa[1, j]);
  }

  vector sai_l2(int j, vector t, real tau, matrix log_thetaa) {
    return (t - tau) * exp(-log_thetaa[2, j])
           + tau    * exp(-log_thetaa[1, j]);
  }

  real ssalt_lpdf(vector t,
                  vector C_i,
                  real tau,
                  int n1,
                  int nc,
                  real beta1,
                  real beta2,
                  matrix log_thetaa) {

    int n = num_elements(t);
    real loglik = 0;

    vector[n] s11 = sai_l1(1, t, log_thetaa);
    vector[n] s12 = sai_l1(2, t, log_thetaa);
    vector[n] s21 = sai_l2(1, t, tau, log_thetaa);
    vector[n] s22 = sai_l2(2, t, tau, log_thetaa);

    // --- Failures before tau ---
    for (i in 1:n1) {
      if (C_i[i] == 1) {
        loglik += log(beta1)
                  - log_thetaa[1,1]
                  + (beta1 - 1) * log(s11[i])
                  - pow(s11[i], beta1)
                  - pow(s12[i], beta2);
      } else {
        loglik += log(beta2)
                  - log_thetaa[1,2]
                  + (beta2 - 1) * log(s12[i])
                  - pow(s11[i], beta1)
                  - pow(s12[i], beta2);
      }
    }

    // --- Failures after tau ---
    for (i in (n1 + 1):nc) {
      if (C_i[i] == 1) {
        loglik += log(beta1)
                  - log_thetaa[2,1]
                  + (beta1 - 1) * log(s21[i])
                  - pow(s21[i], beta1)
                  - pow(s22[i], beta2);
      } else {
        loglik += log(beta2)
                  - log_thetaa[2,2]
                  + (beta2 - 1) * log(s22[i])
                  - pow(s21[i], beta1)
                  - pow(s22[i], beta2);
      }
    }

    // --- Censored observations ---
    for (i in (nc + 1):n) {
      loglik += -pow(s21[i], beta1)
                - pow(s22[i], beta2);
    }

    return loglik;
  }
}


data {
  int<lower=1> n;
  real<lower=0> tau;
  real<lower=0> tc;
  int<lower=0> n1;
  int<lower=0> nc;
  vector[n] C_i;
  vector[n] t;

  real x1;
  real x2;
  real x0;
  real p;

  // --- gamma hyperparameters ---
  vector[2] gphi11;
  vector[2] gphi21;
  vector[2] gphi31;
  vector[2] gphi12;
  vector[2] gphi22;
  vector[2] gphi32;
}

parameters {

  real<lower=1e-9, upper=100> phi_11;
  real<lower=1e-9, upper=100> phi_21;
  real<lower=1e-9, upper=100> phi_31;

  real<lower=1e-9, upper=100> phi_12;
  real<lower=1e-9, upper=100> phi_22;
  real<lower=1e-9, upper=100> phi_32;

}


transformed parameters {

  real beta1 = phi_31;
  real beta2 = phi_32;

  real a1 = log(phi_11) - log(-log(0.999)) / beta1;
  real b1 = -phi_21;

  real a2 = log(phi_12) - log(-log(0.999)) / beta2;
  real b2 = -phi_22;

  // KEY CHANGE: store log(theta) directly — never compute
  // exp(large number). The likelihood is rewritten to use
  // log_thetaa so exp() only appears inside bounded terms.
  matrix[2,2] log_thetaa;
  log_thetaa[1,1] = a1 + b1 * x1;   // log(theta_11)
  log_thetaa[1,2] = a2 + b2 * x1;   // log(theta_12)
  log_thetaa[2,1] = a1 + b1 * x2;   // log(theta_21)
  log_thetaa[2,2] = a2 + b2 * x2;   // log(theta_22)

  // Recover thetaa on original scale for output only
  // (not used inside likelihood)
  matrix[2,2] thetaa;
  thetaa[1,1] = exp(log_thetaa[1,1]);
  thetaa[1,2] = exp(log_thetaa[1,2]);
  thetaa[2,1] = exp(log_thetaa[2,1]);
  thetaa[2,2] = exp(log_thetaa[2,2]);

  // Newton solve for t_p fully on log scale
  real log_tp = log(10);

  for (k in 1:30) {
    real log_th10 = a1 + b1 * x0;
    real log_th20 = a2 + b2 * x0;

    real s1 = exp(beta1 * (log_tp - log_th10));
    real s2 = exp(beta2 * (log_tp - log_th20));

    real f  = log(1 - p) + s1 + s2;
    real df = beta1 * s1 + beta2 * s2;

    log_tp -= f / df;
  }

  real tp = exp(log_tp);

}


model {

  phi_11 ~ gamma(gphi11[1], gphi11[2]);
  phi_21 ~ gamma(gphi21[1], gphi21[2]);
  phi_31 ~ gamma(gphi31[1], gphi31[2]);

  phi_12 ~ gamma(gphi12[1], gphi12[2]);
  phi_22 ~ gamma(gphi22[1], gphi22[2]);
  phi_32 ~ gamma(gphi32[1], gphi32[2]);

  target += ssalt_lpdf(t | C_i, tau, n1, nc,
                        beta1, beta2, log_thetaa);
}
