##### Required Libraries #####
library(nimbleHMC)
library(dplyr)
library(doParallel)
library(HDInterval)
library(mcmcr)

##### LSRRM Models #####
source('Model_Dist.R')
source('Model_InnerProd.R')
source('Model_Proj.R')

##### Utilities #####
source('Utils.R')

##### Setting (TODO) #####
# Data shall be scaled as 1 -- K in a K-point item.
y = read.csv('toy_dta.csv', header=F) %>% as.matrix()

dta = list(y=y) 

Args = list(
  N = nrow(y), # sample size
  V = 2, # dimension of \xi 
  K = 5, # number of categories
  nchain = 3, # number of chains
  thin = 20, # thinning number
  niter = 60000, # number of iterations
  nburnin = 40000, # number of burn-in
  HPD = .95, # HPD interval
  Model = 'Dist' 
  # Dist: latent distance LSRRM
  # InnerProd: Inner product distance LSRRM
  # Proj: Project distance LSRRM
)

J_idx = NULL
for(i in 1:Args[['N']]){
  J_idx = rbind(J_idx, (1:Args[['N']])[-i])
}

constants = list(
  V=Args[['V']], 
  K=Args[['K']], 
  N=Args[['N']], 
  M=rep(0,2), 
  M_d=rep(0,Args[['V']]), 
  Sigma_d=diag(1,Args[['V']]), 
  J_idx=J_idx,
  zero_K = rep(0, Args[['K']])
) 

inits = list()
for(i in 1:Args[['nchain']]){
  inits[[i]] =
    list(
      theta = cbind(scale(apply(y,1,sum)),scale(apply(y,2,sum))),
      tau = seq(-1.5, 1.5, length.out = Args[['K']]-1 ), 
      Sds = 1,
      Ustar = diag(1,2),
      xi = matrix(0,  Args[['N']],  Args[['V']]),
      log_lambda = c(-5, .5),
      U = diag(1,2),
      omega = .1,
      delta = 0
    )
}

##### Compiling #####
my.cluster <- parallel::makeCluster(Args[['nchain']], type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

if(Args[['Model']]=='Dist'){
  chain_output = foreach(c = 1:Args[['nchain']], .combine = list, .multicombine=TRUE) %dopar% {
    
    mod = LSRRM_Dist(data = dta, 
                     constants = constants, inits = inits, 
                     thin = Args[['thin']], niter = Args[['niter']], nburnin = Args[['nburnin']],
                     cal_dist = F) 
    
    out = list(mod)
    return(out)
  }
}else if(Args[['Model']]=='InnerProd'){
  chain_output = foreach(c = 1:Args[['nchain']], .combine = list, .multicombine=TRUE) %dopar% {
    
    mod = LSRRM_InnerProd(data = dta, 
                          constants = constants, inits = inits, 
                          thin = Args[['thin']], niter = Args[['niter']], nburnin = Args[['nburnin']],
                          cal_dist = F) 
    
    out = list(mod)
    return(out)
  }
}else if(Args[['Model']]=='Proj'){
  chain_output = foreach(c = 1:Args[['nchain']], .combine = list, .multicombine=TRUE) %dopar% {
    
    mod = LSRRM_Proj(data = dta, 
                     constants = constants, inits = inits, 
                     thin = Args[['thin']], niter = Args[['niter']], nburnin = Args[['nburnin']],
                     cal_dist = F) 
    
    out = list(mod)
    return(out)
  }
}

stopCluster(my.cluster)

##### Estimating #####
Est = NULL
for(c in 1:Args[['nchain']]){
  Est = rbind(Est, chain_output[[c]][[1]]$summary[,'Mean'])
}
Est = apply(Est, 2, mean)

n_var_est = names(Est)

sel_Tau = grep('tau', n_var_est)
sel_Sigma = grep('Sigma', n_var_est)
sel_Theta = grep('theta', n_var_est)
sel_Omega = grep('omega', n_var_est)
sel_Lambda = grep('lambda', n_var_est)
sel_Xi = grep('xi', n_var_est)

Est_Tau = Est[sel_Tau]
Est_Sigma = Est[sel_Sigma] %>% matrix(2,2)
Est_Theta = Est[sel_Theta] %>% matrix(Args[['N']], 2)
Est_Omega = Est[sel_Omega]
Est_Lambda = Est[sel_Lambda][(Est_Omega>0.5)+1]

##### Procrustes matching #####
Est_Xi = ProcrustesMatching(chain_output, Args=Args)
names(Est_Xi) = n_var_est[sel_Xi]

##### HPD #####
post_samples = NULL
for(c in 1:Args[['nchain']]){
  post_samples = rbind(post_samples, as.matrix(chain_output[[c]][[1]]$samples))
}
hdi_Tau = post_samples[,sel_Tau] %>% hdi(Args[['HPD']])
hdi_Sigma = post_samples[,sel_Sigma] %>% hdi(Args[['HPD']])
hdi_Theta = post_samples[,sel_Theta] %>% hdi(Args[['HPD']])
hdi_Omega = post_samples[,sel_Omega] %>% hdi(Args[['HPD']])
hdi_Lambda = post_samples[,sel_Lambda][,(Est_Omega>0.5)+1] %>% hdi(Args[['HPD']])

##### R-Square #####
rhat_df = rhat(as.mcmc(post_samples), 'parameter', as_df = T)

##### DIC value: LogL: -2logLikelihood; DIC: DIC #####
DIC = DIC_LSRRM(y, nchain=Args[['nchain']], K=Args[['K']], v=Args[['V']], Dist=Args[['Model']])

##### Sender Fit & Receiver Fit #####
Rep_Tau = post_samples[,sel_Tau] %>% array(dim=c(nrow(post_samples), Args[['K']]-1))
Rep_Theta = post_samples[,sel_Theta] %>% array(dim=c(nrow(post_samples), Args[['N']], 2))
Rep_Omega = post_samples[,sel_Omega] 
Rep_Lambda = post_samples[,sel_Lambda] %>% array(dim=c(nrow(post_samples), 2))
Rep_Xi = post_samples[,sel_Xi] %>% array(dim=c(nrow(post_samples), Args[['N']], Args[['V']]))

IN_s = IN_r = NULL
IN_s_p = IN_r_p = NULL
for(i in 1:Args[['N']]){
  # preference
  FIT_s = Fit_s(y[i,-i], Est_Theta[i,1], Est_Theta[-i,2], Est_Tau, Est_Lambda, Est_Xi[i,], Est_Xi[-i,], Mod=Args[['Model']])
  FIT_r = Fit_r(y[-i,i], Est_Theta[-i,1], Est_Theta[i,2], Est_Tau, Est_Lambda, Est_Xi[-i,], Est_Xi[i,], Mod=Args[['Model']])
  tmp_IN_s_V = tmp_IN_r_V = 0
  FIT_rep_s_all = vector('numeric', nrow(post_samples))
  FIT_rep_r_all = vector('numeric', nrow(post_samples))
  for(reps in 1:nrow(post_samples)){
    y_rep_s = data_gen_fit(Rep_Theta[reps,i,1], Rep_Theta[reps,-i,2], Rep_Tau[reps,], Rep_Lambda[reps, (Rep_Omega[reps]>0.5)+1], Rep_Xi[reps,i,], Rep_Xi[reps,-i,], 's', Mod=Args[['Model']])[1,]
    FIT_rep_s = Fit_s(y_rep_s, Rep_Theta[reps,i,1], Rep_Theta[reps,-i,2], Rep_Tau[reps,], Rep_Lambda[reps, (Rep_Omega[reps]>0.5)+1], Rep_Xi[reps,i,], Rep_Xi[reps,-i,], Mod=Args[['Model']])
    tmp_IN_s_V = tmp_IN_s_V + (FIT_rep_s > FIT_s)
    FIT_rep_s_all[reps] = FIT_rep_s
    
    y_rep_r = data_gen_fit(Rep_Theta[reps,-i,1], Rep_Theta[reps,i,2], Rep_Tau[reps,], Rep_Lambda[reps, (Rep_Omega[reps]>0.5)+1], Rep_Xi[reps,-i,], Rep_Xi[reps,i,], 'r', Mod=Args[['Model']])[,1]
    FIT_rep_r = Fit_r(y_rep_r,Rep_Theta[reps,-i,1], Rep_Theta[reps,i,2], Rep_Tau[reps,], Rep_Lambda[reps, (Rep_Omega[reps]>0.5)+1], Rep_Xi[reps,-i,], Rep_Xi[reps,i,], Mod=Args[['Model']])
    tmp_IN_r_V = tmp_IN_r_V + (FIT_rep_r > FIT_r)
    FIT_rep_r_all[reps] = FIT_rep_r
  }
  
  IN_s_p = c(IN_s_p, tmp_IN_s_V / nrow(post_samples))
  IN_r_p = c(IN_r_p, tmp_IN_r_V / nrow(post_samples))
  
  IN_s = c(IN_s, FIT_s)
  IN_r = c(IN_r, FIT_r)
}
# Sender Fit, p_value | Receiver Fit, p_value
cbind(SenderFit=IN_s, SenderFit_P=IN_s_p, Receiver_Fit=IN_r, Receiver_Fit_p=IN_r_p)
