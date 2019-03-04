/*
 Stan model created for Magda
 */

data{
  real alpha_mu;
  real alpha_sigma;
  real beta_mu;
  real beta_sigma;
  int N;
  int n_grps;
  int grp[N];
  vector [N] y;
}

parameters{
  vector[n_grps] log_mu;
  vector<lower = 0>[n_grps] log_sigma;
}

model{
  log_mu ~ normal(alpha_mu, beta_mu);
  log_sigma ~ normal(alpha_sigma, beta_sigma);
  y ~ lognormal(log_mu[grp], log_sigma[grp]);
}

generated quantities{
  vector[n_grps] mu; 
  vector[n_grps] sigma; 
  for(n in 1:n_grps){
	mu[n] = exp(log_mu[n] + 0.5 * log_sigma[n]^2);
	sigma[n] = sqrt(mu[n]^2 * (exp(log_sigma[n]^2) - 1));
  }
  
}

