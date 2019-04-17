
      model{
      # Priors
      
      mu_log ~ dnorm(0, 0.001)
      sigma_log ~ dunif(0,60)
      
      mu_y= exp(mu_log+sigma_log^2/2)
      sigma_y = sqrt(exp(2*mu_log+sigma_log^2)*(exp(sigma_log^2)-1))
      tau_y = 1/(sigma_y^2)
      
      # Likelihood
      for(i in 1:n){
      mu[i] = mu_group[y.POP.group[i]] 
      sigma[i] = sigma_group[y.POP.group[i]]
      tau [i] = 1/sigma [i]^2
      y.SLA[i] ~ dnorm(mu[i], tau [i]) }

      for (j in 1:y.POP.n){ # set loop to cycle through POPS
      mu_group[j] ~ dlnorm(mu_y, tau_y) # To make pops TOTALLY
      sigma_group[j] ~ dunif(0, 60) 
      }}
      
