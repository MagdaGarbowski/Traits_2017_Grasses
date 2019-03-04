/*
 Stan model created for Magda
 */

data{
  real mean_mu;
  real mean_sigma;
  real var_mu;
  real var_sigma;
  real var_H_num;
  int N;
  int n_grps;
  int grp[N];
  int H_num[N];
  vector[N] y;
}

parameters{
  vector[3] beta_H_num;
  vector[n_grps] log_mu;
  vector<lower = 0>[n_grps] log_sigma;
}

model{
  // Priors
  beta_H_num ~ normal(0, var_H_num);
  log_mu ~ normal(mean_mu, mean_sigma);
  log_sigma ~ normal(var_mu, var_sigma);

  // Likelihood
  y ~ lognormal(log_mu[grp] + beta_H_num[H_num], log_sigma[grp]);
}

generated quantities{
  vector[n_grps] mu[max(H_num)]; 
  vector[n_grps] sigma[max(H_num)]; 
  for(t in 1:max(H_num))
	for(n in 1:n_grps){
	  mu[t,n] = exp(log_mu[n] + beta_H_num[t] + 0.5 * log_sigma[n]^2);
	  sigma[t,n] = sqrt(mu[t,n]^2 * (exp(log_sigma[n]^2) - 1));
	}
  
}

