library(psych)
library(pscl)
library(dplyr)
library(mvtnorm)
library(foreach)
library(doParallel)
library(ranger)
library(tidyverse)

seqlast = function(from, to, by){
  vec = do.call(seq, list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

SNIRT_UND_RSM = function(y, J, nchain=2, Update=F, mod=NULL, epoch=300, nbatch=min(c(nrow(y),30)), iters=(epoch*N/nbatch)+1, burnin=round((iters-1)/3), thin=round((iters-burnin-1)/1000)){
  
  # MCMC
  sb = function(x,a,b){
    new_b = rgamma(1,2.5+a,x*b/(x+b))
    return(new_b)
  }
  
  ss3 = function(x,a,b){
    I = length(x)
    s = pscl::rigamma(1,a+I/2,b+sum(x^2)/2)
    return(s)
  }
  
  ptau = function(y,tp,tf,tau,St){
    
    N = length(tp)
    J = length(tau) +1
    
    
    prob = array(1, dim=c(N,N,J))
    m_tp = matrix(tp,N,N, byrow=F) 
    m_tf = matrix(tf,N,N, byrow=T)
    m_t = array(m_tp + m_tf, dim=c(N,N,J-1))
    
    m_I = array(rep(1:(J-1),each=N^2), dim=c(N,N,J-1))
    mc_tau = array(rep(cumsum(tau),each = N^2),dim=c(N,N,J-1) )
    prob[,,2:J] = exp((m_I *m_t - mc_tau) ) 
    c_Prob = array(rowSums(prob, dims = 2), dim=c(N,N,J))
    
    m_I = array(rep(1:J,each=N^2), dim=c(N,N,J))  
    m_y = array(y, dim=c(N,N,J))
    Prob = (m_y==m_I)*prob/c_Prob
    sum_Prob = rowSums(Prob, dims=2)
    diag(sum_Prob) = 1
    
    (-sum(tau^2)/(2*St))+sum(log(sum_Prob))
  }
  
  ptheta = function(y,TP,tf,tau,Stp){
    
    N = length(tf) 
    J = length(tau) +1
    
    ## Tp
    tpf = (TP + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp(m_tpf - mc_tau)
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = matrix(y, nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob_p = apply(Prob * m_y,1,sum)
    
    (-2*0.5*(sum(TP^2)/Stp)+2*sum(log(sum_Prob_p)))
  }
  
  pSTT = function(STT, T1, a=2.5, b){
    N = length(T1)
    prior_tmp = log((1/STT)^(a+1)) + (-b/STT)  
    likeli_tmp = log((1/STT)^(N/2)) + (-sum(T1^2)/(2*STT)) 
    
    prior_tmp + likeli_tmp
  } 
  
  # Estimate
  eap = function(x,node){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,node,c]/nchain 
    }
    out = apply(tmp, 1, mean)
    return(out)
  } 
  
  eap1d = function(x,node){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,node,c]/nchain 
    }
    out = mean(tmp)
    return(out)
  } 
  
  MCMC = function(x){
    tmp = 0
    for(c in 1:nchain){
      tmp = tmp + x[,,c]/nchain 
    }
    return(tmp)
  } 
  
  RR_s = function(x){
    len = dim(x)[1]
    nchain = dim(x)[2]
    B = NULL
    W = NULL
    B_tmp = NULL
    W_tmp = NULL
    for(c in 1:nchain){
      B_tmp = c(B_tmp, mean(x[,c]))
      W_tmp = c(W_tmp, var(x[,c]))
    }
    B = c(B, len*var(B_tmp))
    W = c(W, mean(W_tmp))
    
    Var = (len-1)*W/len + B/len
    sqrt(Var/W)
  }
  
  RR = function(x){
    items = dim(x)[1]
    len = dim(x)[2]
    nchain = dim(x)[3]
    B = NULL
    W = NULL
    for(i in 1:items){
      B_tmp = NULL
      W_tmp = NULL
      for(c in 1:nchain){
        B_tmp = c(B_tmp, mean(x[i,,c]))
        W_tmp = c(W_tmp, var(x[i,,c]))
      }
      B = c(B, len*var(B_tmp))
      W = c(W, mean(W_tmp))
    }
    Var = (len-1)*W/len + B/len
    sqrt(Var/W)
  }
  
  r = 0
  N = nrow(y)
  K = J
  Batch = seqlast(1,N+1,by=nbatch)
  BatchLen = length(Batch)-1
  
  ## Initial value
  {
    ### Chain
    tp = array(scale(apply(y,1,sum)),dim = c(N, iters, nchain))
    tf = scale(apply(y,2,sum))
    tau = array(seq(-2,2,length.out=J-1),dim = c(J-1, iters, nchain))
    Stp = array(1,dim = c(1, iters, nchain))
    St = 0.5
  }
  
  {
    Jtau = rep(0.1,J-2)
    Jtp = rep(0.5,N)
    Jstp = 0.5
  }
  
  # Parallel backend
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(min(n.cores,nchain), type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  ## Estimation
  tmp = foreach(c = 1:nchain, .combine = list, .multicombine=TRUE) %dopar%{
    iter = 1
    for(e in 1:epoch){
      for(b in 1:BatchLen){
        posi = sample(1:N, replace = F)
        batch_seq = Batch[b]:(Batch[b+1]-1) 
        batch = posi[batch_seq]
  
        ### tau
        {
          proposal = mvtnorm::rmvnorm(1,tau[1:(J-2),iter,c],diag(Jtau,J-2))
          proposal = sort(c(proposal,-sum(proposal)))
          R = exp(ptau(y[batch,batch],tp[batch,iter,c],tf[batch],proposal,St)
                 -ptau(y[batch,batch],tp[batch,iter,c],tf[batch],tau[,iter,c],St) )
          if(is.na(R)==T){R=0}
          accept = rbinom(1,1,min(1,R))
          tau[,iter+1,c] = ifelse(rep(accept,K-1),proposal,tau[,iter,c])
          
          Jtau = apply(tau[,1:(iter+1),c],1,var)[1:(J-2)] + rep(0.0001,J-2)
        }
        
        ### tau
        {
          St = ss3(tau[,iter+1,c],3,10)
        }
        
        bstp = sb(Stp[,iter,c], 2.5, 2)
        
        ### Stp 
        {
          ### Stp
          proposal = rlnorm(1,log(Stp[,iter,c]),Jstp)
          R = exp(log(proposal)+pSTT(proposal,         tp[,iter,c],2.5,bstp)
                  -log(Stp[,iter,c])-pSTT(Stp[,iter,c],tp[,iter,c],2.5,bstp) )
          R[is.na(R)] = 0
          accept = rbinom(1,1,min(1,R))
          Stp[,iter+1,c] = ifelse(accept,proposal,Stp[,iter,c])
          
          Jstp = sd(log(Stp[,1:(iter+1),c]))+0.00001
        }
        
        for(n in 1:N){
          n_theta = unique(c(n, batch))[-1]
          {
            ### theta
            proposal = rnorm(1,tp[n,iter,c],1.5)
            R = exp(ptheta(y[n,n_theta],proposal,tf[n_theta],tau[,iter+1,c],Stp[,iter+1,c])-
                ptheta(y[n,n_theta],tp[n,iter,c],tf[n_theta],tau[,iter+1,c],Stp[,iter+1,c]))
            if(is.na(R)==T){R=0}
            accept = rbinom(1,1,min(1,R))
            tp[n,iter+1,c] = ifelse(accept,proposal,tp[n,iter,c])
            tf[n] = tp[n,iter+1,c]
          }
        }
        
        iter = iter + 1
      }
    }
    out = list(tau = tau[,,c],
               tp = tp[,,c],
               Stp = Stp[,,c])
    return(out)
    
  }
  
  stopCluster(my.cluster)
  
  for(i in 1:nchain){
    tp[,,i] = tmp[[i]]$tp
    tau[,,i] = tmp[[i]]$tau
    Stp[,,i] = tmp[[i]]$Stp
  }
  
  # EAP
  node = seq(burnin+thin,iters,by=thin)
  
  ID = row.names(y)
  
  tp_out = eap(tp, node); names(tp_out) = ID
  tau_out = eap(tau, node) 
  Stp_out = eap1d(Stp, node)
  
  # tp_MCMC = MCMC(tp)
  # tf_MCMC = MCMC(tf)
  # tau_MCMC = MCMC(tau)
  
  tp_RR = RR(tp[,node,])
  tau_RR = RR(tau[,node,])
  Stp_RR = RR_s(Stp[,node,])
  
  out = list(
    Estimate = list(
      Theta = tp_out,
      Tau = tau_out,
      Stp = Stp_out
    ),
    MCMC = list(
      Theta = tp,
      Tau = tau,
      Stp = Stp
    ),
    RSquare = list(
      Theta = tp_RR,
      Tau = tau_RR,
      Stp = Stp_RR
    ),
    Node = node 
  )
  return(out)
}

