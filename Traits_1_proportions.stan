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
  vector[n_grps] logit_mu;
  vector<lower = 0>[n_grps] kappa;
}

model{
  // Priors
  beta_H_num ~ normal(0, var_H_num);
  logit_mu ~ normal(mean_mu, mean_sigma);  // change these to logistic as well? 
  kappa ~ normal(var_mu, var_sigma);

  // Likelihood
  y ~ beta(inv_logit(logit_mu[grp] + beta_H_num[H_num]).*kappa[grp], (1-inv_logit(logit_mu[grp] + beta_H_num[H_num])).*kappa[grp]); // Look up element-wise multiplication in STAN
}

generated quantities{
  vector[n_grps] mu[max(H_num)]; 
  vector[n_grps] sigma[max(H_num)]; 
  vector[n_grps] cv[max(H_num)]; 
  for(t in 1:max(H_num))
	for(n in 1:n_grps){
	  mu[t,n] = inv_logit(logit_mu[n]+ beta_H_num[t]); // To back transform use inv_logit (logistic_mu[n])
	  sigma[t,n] = (mu[t,n] * (1-(mu[t,n])))/(1+kappa[n]) ; // this will be a formula with mu and kappa. 
	  cv[t,n] = sigma[t,n] / mu[t,n];
	}
  
}

