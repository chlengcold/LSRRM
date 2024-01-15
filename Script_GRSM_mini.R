library(psych)
library(pscl)
library(dplyr)
library(mvtnorm)
library(foreach)
library(doParallel)
library(ranger)
library(tidyverse)

SNIRT_UND_GRSM = function(y, J, nchain=2, iters=4500, burnin=1500, thin=3, doc1){
  
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
  
  ss2 = function(x,a,b){
    I = length(x)
    s = pscl::rigamma(1,a+I,b+sum(x))
    return(s)
  }
  
  ptau = function(y,a,tp,tf,tau,St){
    
    N = length(tp)
    J = length(tau) +1
    
    prob = array(1, dim=c(N,N,J))
    m_tp = matrix(tp,N,N, byrow=F) 
    m_tf = matrix(tf,N,N, byrow=T)
    m_t = array(a*m_tp+m_tf, dim=c(N,N,J-1))
    
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
  
  ptheta = function(y,TP,a,tp,tf,tau){
    
    N = length(tf) 
    J = length(tau) +1
    
    ## Tp
    tpf = (a*TP + tf)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf - mc_tau))
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = matrix(y[1,-1], nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob_p = apply(Prob * m_y,1,sum)
    
    ## Tp
    tpf = (a*tp + TP)
    m_I = matrix(0:(J-1), N, J, byrow=T)
    m_tpf = matrix(tpf, N, J, byrow=F)*m_I
    c_tau = cumsum(tau)
    mc_tau = matrix(c_tau, N, J-1, byrow=T)
    mc_tau = cbind(0, mc_tau)
    prob = exp((m_tpf - mc_tau))
    c_Prob = matrix(apply(prob, 1, sum), N,J, byrow=F)
    Prob = prob / c_Prob
    m_y = matrix(y[1,-1], nrow=N, ncol=J, byrow=F) == matrix(1:J, nrow=N, ncol=J, byrow=T)
    sum_Prob_p = apply(Prob * m_y,1,sum)
    
    (-(sum(TP^2))+sum(log(sum_Prob_p)))
  }

  pa = function(y,a,tp,tf,tau,Sa){

    N = length(tp)
    J = length(tau) +1
    
    prob = array(1, dim=c(N,N,J))
    Prob = array(0, dim=c(N,N,J))
    c_Prob = matrix(1,N,N)
    sum_Prob = matrix(0,N,N)
    m_tp = matrix(tp,N,N, byrow=F) 
    m_tf = matrix(tf,N,N, byrow=T)
    m_t = (m_tp + m_tf)
    for(j in 2:J){
      mc_tau = matrix(sum(tau[1:(j-1)]),N,N)
      prob[,,j] = exp( a * ((j-1)*m_t - mc_tau )) 
      c_Prob = c_Prob + prob[,,j] 
    }
    for(j in 1:J){
      Prob[,,j] = (y==j)*(prob[,,j] / c_Prob) 
      tmp = Prob[,,j]
      sum_Prob = sum_Prob + tmp
    }
    diag(sum_Prob) = 1

    -a/Sa+sum(log(sum_Prob))
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
  
  RR = function(x){
    n = length(x)/2
    x1 = x[1:n]
    x2 = x[(n+1):(2*n)]
    B = n*var(c(mean(x1),mean(x2)))
    W = (var(x1)+var(x2))/2
    Var = (n-1)*W/n+B/n
    
    sqrt(Var/W)
  }
  
  
  r = 0
  N = nrow(y)
  K = J
  
  ## Initial value
  {
    ### Chain
    a = array(1,dim = c(1, iters, nchain))
    tp = array(0,dim = c(N, iters, nchain))
    tf = array(0,dim = c(N, iters, nchain))
    tau = array(seq(-2,2,length.out=J-1),dim = c(J-1, iters, nchain))
    St = 0.5
    Sa = 0.5
  }
  
  time1 = Sys.time()
  
  # Parallel backend
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()
  foreach::getDoParWorkers()
  
  ## Estimation
  tmp = foreach(c = 1:nchain, .combine = list) %dopar%{
    for(iter in 1:(iters-1)){
  
      ### tau
      {
        proposal = mvtnorm::rmvnorm(1,tau[1:(J-2),iter,c],diag(0.1,J-2))
        proposal = (c(proposal,-sum(proposal)))
        R = exp(ptau(y,a[,iter,c], tp[,iter,c],tf[,iter,c],proposal,St)
                -ptau(y,a[,iter,c], tp[,iter,c],tf[,iter,c],tau[,iter,c],St) )
        if(is.na(R)==T){R=0}
        accept = rbinom(1,1,min(1,R))
        tau[,iter+1,c] = ifelse(rep(accept,K-1),proposal,tau[,iter,c])
      }
      
      ### a
      proposal = rlnorm(1,log(a[,iter,c]),0.01)
      R = exp(log(proposal)+pa(y,proposal,tp[,iter,c],tf[,iter,c],tau[,iter+1,c],Sa)
              -log(a[,iter,c])-pa(y,a[,iter,c],tp[,iter,c],tf[,iter,c],tau[,iter+1,c],Sa) )
      R[is.na(R)] = 0
      accept = rbinom(1,1,min(1,R))
      a[,iter+1,c] = ifelse(accept,proposal,a[,iter,c])
      
      ### var
      {
        St = ss3(tau[,iter+1,c],3,10)
        Sa = ss2(a[,iter+1,c],3,8)
      }
      
      n_steps = sample(1:N, N, replace = F)
      tp_tmp = tp[,iter,c]
      tf_tmp = tf[,iter,c]
      for(n in n_steps){
        {
          ### theta
          proposal = rnorm(1,tp_tmp[n],1.5)
          R = exp(ptheta(y,proposal,a[,iter+1,c],tp_tmp,tf_tmp,tau[,iter+1,c],n)-
                 ptheta(y,tp_tmp[n],a[,iter+1,c],tp_tmp,tf_tmp,tau[,iter+1,c],n))
          if(is.na(R)==T){R=0}
          accept = rbinom(1,1,min(1,R))
          tp_tmp[n] = ifelse(accept,proposal,tp_tmp[n])
          tf_tmp[n] = ifelse(accept,proposal,tf_tmp[n])
        }
      }
      tf[,iter+1,c] = tf_tmp
      tp[,iter+1,c] = tp_tmp
      
    }
    out = list(tau = tau[,,c],
               a = a[,,c],
               tp = tp[,,c],
               tf = tf[,,c])
    return(out)
    
  }
  
  stopCluster(my.cluster)
  
  for(i in 1:nchain){
    tp[,,i] = tmp[[i]]$tp
    tf[,,i] = tmp[[i]]$tf
    a[,,i] = tmp[[i]]$a
    tau[,,i] = tmp[[i]]$tau
  }
  
  # EAP
  node = seq(burnin+thin,iters,by=thin)
  
  ID = row.names(y)
  
  tp_out = eap(tp, node); names(tp_out) = ID
  tf_out = eap(tf, node); names(tf_out) = ID 
  a_out = eap1d(a, node); names(a_out) = ID 
  tau_out = eap(tau, node) 
  
  tp_RR = apply(cbind(tp[,,1],tp[,,2]),1,RR)
  tf_RR = apply(cbind(tf[,,1],tf[,,2]),1,RR)
  a_RR = RR(cbind(a[,,1],a[,,2]))
  tau_RR = apply(cbind(tau[,,1],tau[,,2]),1,RR)
  
  write.table(round(tp_out,4),paste0(doc1,"\\tp.txt"),col.names = F,row.names = T,sep="\t")
  write.table(round(tf_out,4),paste0(doc1,"\\tf.txt"),col.names = F,row.names = T,sep="\t")
  write.table(round(a_out,4),paste0(doc1,"\\a.txt"),col.names = F,row.names = T,sep="\t")
  write.table(round(tau_out,4),paste0(doc1,"\\tau.txt"),col.names = F,row.names = F,sep="\t")
  
  for(c in 1:nchain){
    write.table(round(a[,,c],4),paste0(doc1,"\\a_MCMC_",c,".txt"),col.names = F,row.names = F,sep="\t")
    write.table(round(tp[,,c],4),paste0(doc1,"\\tp_MCMC_",c,".txt"),col.names = F,row.names = F,sep="\t")
    write.table(round(tf[,,c],4),paste0(doc1,"\\tf_MCMC_",c,".txt"),col.names = F,row.names = F,sep="\t")
    write.table(round(tau[,,c],4),paste0(doc1,"\\tau_MCMC_",c,".txt"),col.names = F,row.names = F,sep="\t")
  }

  out = list(
    Estimate = list(
      ThetaP = tp_out,
      ThetaF = tf_out,
      A = a_out,
      Tau = tau_out
    ),
    MCMC = list(
      ThetaP = tp,
      ThetaF = tf,
      A = a,
      Tau = tau
    ),
    RSquare = list(
      ThetaP = tp_RR,
      ThetaF = tf_RR,
      A = a_RR,
      Tau = tau_RR
    )
  )
  return(out)
}

