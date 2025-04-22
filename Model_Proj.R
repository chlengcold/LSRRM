# Defining model
LSRRM_Proj = function(data, constants, inits, thin = 40, niter = 60000, nburnin = 20000, cal_dist=F) {
  library(nimbleHMC)
  
  code = nimbleCode({ 
    # Model
    for(i in 1:N) {   
      dist[i, i] <- 0 
      for(j in J_idx[i,]) {	 
        dist[i, j] <- sum(xi[i,1:V] * xi[j,1:V]) / sqrt(sum(xi[j,1:V]^2))  

        l[ i, j, 1] <- 1		
        for(k in 2:K){
          l[i, j, k] <- exp((k-1)*(theta[i,1] + theta[j,2] + lambda[delta+1] * dist[i,j]) - sum(tau[1:(k-1)]) ) 
        }
        
        cl[i,j] <- sum(l[i, j, 1:K]) 
        
        prob[i, j, 1:K] <- l[i, j, 1:K]/cl[i,j]
        
        y[i , j] ~ dcat(prob[i,j,1:K]) 
      } 
    }                              
    
    # Priors  
    log_lambda[1] ~ dnorm(-5,1)
    log_lambda[2] ~ dnorm(.5,1)
    lambda[1:2] <- exp(log_lambda[1:2])
    omega ~ dbeta(1,1)
    delta ~ dbern(omega) 
    
    for(k in 1:(K-1)){
      tau_raw[k] ~ dnorm(0, var=4)
    }
    tau[1:(K-1)] <- tau_raw[1:(K-1)] - mean(tau_raw[1:(K-1)])
    
    
    Ustar[1:2,1:2] ~ dlkj_corr_cholesky(1,2) # upper-triangular
    U[1:2, 1:2] <- t(Ustar[1:2, 1:2]) %*% Ustar[1:2, 1:2] # transpose to lower times itself
    Sds ~ T(dt(0, tau=1/2.5, 1),0,)
    Sigma[1,1] <- Sds
    Sigma[2,2] <- Sds
    Sigma[1,2] <- U[2,1] * Sds
    Sigma[2,1] <- U[2,1] * Sds
    
    for(i in 1: N){
      theta[i,1:2] ~ dmnorm(M[1:2], cov=Sigma[1:2, 1:2])
      xi_raw[i,1:V] ~ dmnorm(M_d[1:V], cov=Sigma_d[1:V, 1:V])
      
      logL_raw[i, i, 1:K] <- zero_K[1:K]
      for(j in J_idx[i,]){
        for(k in 1:K){
          logL_raw[i, j, k] <- log(prob[i, j, k]) * (y[i, j] == k)
        } # k
      } # j
    }
    
    for(v in 1:V){
      xi[1:N, v] <- xi_raw[1:N, v] - mean(xi_raw[1:N, v])
    }
  })
  Rmodel <- nimbleModel(code,  data = dta, inits = inits, constants = constants, buildDerivs = TRUE)
  
  ## create default MCMC configuration object 
  monitors = c('theta','xi', 'omega','delta', 'lambda', 'tau', 'Sigma', 'dist', 'logL_raw')
  if(cal_dist==F){
    monitors = monitors[-grep('dist', monitors)]
  }
  conf <- configureMCMC(Rmodel, enableWAIC = F, monitors = monitors) 
  
  ## remove default samplers 
  addHMC(conf, target = c('Ustar', 'theta', 'log_lambda', 'omega', 'xi_raw', 'tau_raw'), replace=T) 
  Rmcmc <- buildMCMC(conf)
  
  Cmodel <- compileNimble(Rmodel) 
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel) 
  samples <- runMCMC(Cmcmc, thin = thin, niter = niter, nburnin = nburnin, nchains = 1,
                     summary = TRUE, WAIC = F, samplesAsCodaMCMC=T)
  
  return(samples)
}
