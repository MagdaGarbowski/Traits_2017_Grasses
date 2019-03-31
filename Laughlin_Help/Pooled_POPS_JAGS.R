
model{
# Priors
  mu_log_y ~ dnorm(0, 0.001)
  sigma_log_y ~ dunif(0,60)
  mu_Y <- exp(mu_log_y+sigma_log_y^2/2)
  sigma_Y <- sqrt(exp(2*mu_log_y+sigma_log_y^2)*(exp(sigma_log_y^2)-1))
        # Likelihood
        for(i in 1:n){
        y.SLA[i]~dlnorm(mu_log_y,1/sigma_log_y^2)}
      }
  
